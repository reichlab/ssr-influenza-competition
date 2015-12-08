library(ssr)
library(magrittr)
library(plyr)
library(dplyr)

## all functions in ssr.R
##   leading X means test written,
##   leading I means implicitly tested by another test
##   leading S means simple enough that no test is required
##   no leading character means test needs to be written still
# S "assemble_prediction_examples"                       
# X "assemble_training_examples"                         
# X "compute_kernel_values"                              
# X "compute_lagged_obs_vecs"                            
# X "compute_normalized_log_weights"                     
# S "create_ssr_control"                                 
# S "create_ssr_control_default"                         
#   "discrete_kernel"                                    
#   "est_ssr_params_stepwise_crossval"                   
#   "est_ssr_params_stepwise_crossval_one_potential_step"
# S "get_default_kernel_fns"                             
#   "get_dist_predictions_one_week"                      
# S "get_inds_smallest_k"                                
# X "get_kernel_fn_init_params"                          
#   "get_prediction_season_week"                         
#   "get_pt_predictions_one_week"                        
# S "periodic_kernel"                                    
# S "squared_exp_kernel"                                 
#   "ssr"                                                
#   "ssr_crossval_estimate_parameter_loss"               
#   "ssr_predict"                                        
#   "ssr_predict_dengue_one_week"                        
#   "ssr_predict_given_lagged_obs"                       
# X "unvectorize_theta"                                  
# I "unvectorize_theta_one_kernel_fn"                    
#   "validate_ssr_control"                               
# S "vectorize_theta"



context("ssr prediction functions")

test_that("get_kernel_fn_init_params works", {
    expected <- list(bw=1, period=1)
    
    actual <- get_kernel_fn_init_params("b",
        ssr_control=list(theta_est=list(a="bw", b=c("bw", "period"))))
    
    expect_identical(actual, expected)
})

test_that("unvectorize_theta works, no fixed params", {
    theta <- list(a_lag1=list(bw=1),
        a_lag5=list(bw=15),
        b_lag0=list(bw=3,
            period=12),
        b_lag2=list(bw=2,
            period=987))
    
    lags <- list(a=c(1, 5),
        b=c(0, 2))
    
    expected <- theta
    
    actual <- unvectorize_theta(vectorize_theta(theta),
        lags,
        ssr_control=list(theta_est=list(a="bw", b=c("bw", "period"))),
        add_fixed_params=FALSE)
    
    expect_identical(actual, expected)
})

test_that("unvectorize_theta works, fixed params not included", {
    theta <- list(a_lag1=list(bw=1),
        a_lag5=list(bw=15),
        b_lag0=list(bw=3,
            period=12),
        b_lag2=list(bw=2,
            period=987))
    
    lags <- list(a=c(1, 5),
        b=c(0, 2))
    
    expected <- list(a_lag1=list(bw=1),
        a_lag5=list(bw=15),
        b_lag0=list(bw=3),
        b_lag2=list(bw=2))
    
    actual <- unvectorize_theta(vectorize_theta(theta),
        lags,
        ssr_control=list(theta_est=list(a="bw", b=c("bw")),
            theta_fixed=list(b=list(period=15))),
        add_fixed_params=FALSE)
    
    expect_identical(actual, expected)
})

test_that("unvectorize_theta works, fixed params included", {
    theta <- list(a_lag1=list(bw=1),
        a_lag5=list(bw=15),
        b_lag0=list(bw=3,
            period=12),
        b_lag2=list(bw=2,
            period=987))
    
    lags <- list(a=c(1, 5),
        b=c(0, 2))
    
    expected <- list(a_lag1=list(bw=1),
        a_lag5=list(bw=15),
        b_lag0=list(bw=3, period=15),
        b_lag2=list(bw=2, period=15))
    
    actual <- unvectorize_theta(vectorize_theta(theta),
        lags,
        ssr_control=list(theta_est=list(a="bw", b=c("bw")),
            theta_fixed=list(b=list(period=15))),
        add_fixed_params=TRUE)
    
    expect_identical(actual, expected)
})


test_that("compute_normalized_log_weights works", {
    init_val <- c(0, -100, -50, -78, -99, -0.2)
    
    result <- compute_normalized_log_weights(init_val)
    
    expect_equal(sum(exp(result)), 1)
    expect_true(all.equal(init_val - result,
        rep(init_val[1] - result[1], length(init_val))))
})

test_that("assemble_training_examples works", {
    test_data <- data.frame(a = 1:20, b = rnorm(20))
    lags <- list(a = c(1, 5), b = c(0, 2))
    leading_rows_to_drop <- 6

    actual <- assemble_training_examples(test_data,
        lags,
        y_names="b",
        leading_rows_to_drop=leading_rows_to_drop,
        additional_rows_to_drop=c(14, 15, 18),
        prediction_horizon = 2,
        drop_trailing_rows=TRUE)
    
    expected_lagged <- compute_lagged_obs_vecs(test_data,
        lags,
        c(seq_len(leading_rows_to_drop), c(14, 15, 18, 19, 20)))
    expected_lead <- test_data[1:20 + 2, "b", drop=FALSE]
    expected_lead <- expected_lead[
        -c(seq_len(leading_rows_to_drop), c(14, 15, 18, 19, 20)), , drop=FALSE]
    expected <- list(lagged_obs=expected_lagged,
        lead_obs=expected_lead
    )
    
    expect_identical(actual, expected)
})

test_that("compute_lagged_obs_vecs works -- all variables used", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
    lags <- list(a = c(1, 5), b = c(0, 2))
    leading_rows_to_drop <- 6
    trailing_rows_to_drop <- 1
    
    expected <- test_data
    expected <- expected %>% mutate(a_lag1=lag(a, 1),
        a_lag5=lag(a, 5),
        b_lag0=b,
        b_lag2=lag(b, 2))
    expected <- expected[seq(from=7, to=9), 3:6]
    
    actual <- compute_lagged_obs_vecs(test_data,
        lags,
        c(seq_len(leading_rows_to_drop),
            seq(from=nrow(test_data) - trailing_rows_to_drop + 1,
                to=nrow(test_data))))
    
    expect_identical(actual, expected)
})

test_that("compute_lagged_obs_vecs works -- one variable not used", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
    lags <- list(b = c(0, 2))
    leading_rows_to_drop <- 6
    trailing_rows_to_drop <- 0
    
    expected <- test_data
    expected <- expected %>% mutate(b_lag0=b,
        b_lag2=lag(b, 2))
    expected <- expected[seq(from=7, to=10), 3:4]
    
    actual <- compute_lagged_obs_vecs(test_data,
        lags,
        seq_len(leading_rows_to_drop))
    
    expect_identical(actual, expected)
})


test_that("compute_kernel_values works", {
    test_data <- data.frame(a = 1:10, b = rnorm(10))
    lags <- list(a = c(1), b = c(0, 2, 3))
    leading_rows_to_drop <- 3
    trailing_rows_to_drop <- 1
    
    train_lagged_obs <- compute_lagged_obs_vecs(test_data,
        lags,
        c(seq_len(leading_rows_to_drop),
            seq(from=nrow(test_data) - trailing_rows_to_drop + 1,
                to=nrow(test_data))))
    
    prediction_lagged_obs <- compute_lagged_obs_vecs(test_data,
        lags,
        seq_len(9))
    
    bws <- c(1.3, 1.4, 1.5, 0.8)
    expected <- sapply(seq_len(ncol(train_lagged_obs)), function(ind) {
        squared_exp_kernel(x=train_lagged_obs[, ind],
            center=prediction_lagged_obs[, ind],
            bw=bws[ind],
            log=TRUE)
    })
    expected <- apply(expected, 1, sum)
    
    actual <- compute_kernel_values(train_lagged_obs,
        prediction_lagged_obs,
        lags,
        theta=list(a_lag1=list(bw=1.3),
            b_lag0=list(bw=1.4),
            b_lag2=list(bw=1.5),
            b_lag3=list(bw=0.8)),
        ssr_control=list(kernel_fns=list(a="squared_exp_kernel",
            b="squared_exp_kernel")),
        log = TRUE)
    
    expect_identical(actual, expected)
})

# test_that("ssr_predict with squared exponential kernel works", {
#     test_data <- data.frame(a = 1:10, b = rnorm(10))
#     lags <- list(a = c(1), b = c(0, 2, 3))
#     leading_rows_to_drop <- 3
#     trailing_rows_to_drop <- 1
#     
#     train_lagged_obs <- compute_lagged_obs_vecs(test_data,
#         lags,
#         leading_rows_to_drop,
#         trailing_rows_to_drop)
#     
#     prediction_lagged_obs <- compute_lagged_obs_vecs(test_data,
#         lags,
#         9,
#         0)
#     
#     bws <- c(1.3, 1.4, 1.5, 0.8)
#     expected <- sapply(seq_len(ncol(train_lagged_obs)), function(ind) {
#         squared_exp_kernel(x=train_lagged_obs[, ind],
#             center=prediction_lagged_obs[, ind],
#             bw=bws[ind],
#             log=TRUE)
#     })
#     expected <- apply(expected, 1, sum)
#     
#     actual <- compute_kernel_values(train_lagged_obs,
#         prediction_lagged_obs,
#         lags,
#         theta=list(a_lag1=list(bw=1.3),
#             b_lag0=list(bw=1.4),
#             b_lag2=list(bw=1.5),
#             b_lag3=list(bw=0.8)),
#         control=list(kernel_fns=list(a="squared_exp_kernel",
#             b="squared_exp_kernel")),
#         log = TRUE)
#     
#     expect_identical(actual, expected)
#     ## slow manually computed result that I know is right
#     lag <- 2
#     prediction_horizon <- 3
#     theta <- 2
#     
#     train_data <- 1:12
#     predict_data <- 2:7
#     lagged_train_data <- cbind(1:10, 2:11, 3:12)
#     lagged_predict_data <- cbind(2:5, 3:6, 4:7)
#     
#     dists <- matrix(NA, nrow=nrow(lagged_train_data), ncol=nrow(lagged_predict_data))
#     for(i in seq_len(nrow(lagged_train_data))) {
#         for(j in seq_len(nrow(lagged_predict_data))) {
#             dists[i, j] <- sqrt(sum((lagged_train_data[i, ] - lagged_predict_data[j, ])^2))
#         }
#     }
#     
#     ## cannot get predictions for the first train_lag elements of predict_data:
#     ## not enough preceeding observations to form lagged obs. vectors
#     ## expect NA weights
#     expected_weights <- matrix(NA, nrow=length(train_data), ncol=length(predict_data))
#     
#     ## to predict prediction_step steps ahead from the remaining values of predict_data,
#     for(j in seq(from=train_lag + 1, to=length(predict_data))) {
#         j_lagged <- j - train_lag  # j indexes predict_data, j_lagged indexes lagged_predict_data
#         
#         dists_for_nonzero_weights <- dists[seq(from=1, length=length(train_data) - train_lag - prediction_step), j_lagged]
#         
#         weights_j <- c(rep(0, train_lag), # no weight to the first train_lag elements of train_data -- there was not enough data preceeding them to form lagged obs. vectors
#             rep(0, prediction_step), # no weight to next prediction_step elements of train_data -- the weights that were assigned to the corresponding lagged obs. vectors are projected forward prediction_step steps
#             exp(-theta * dists_for_nonzero_weights / mean(dists_for_nonzero_weights))
#         )
#         
#         weights_j <- weights_j / sum(weights_j)
#         expected_weights[, j] <- weights_j
#     }
# 
#     ssr_fit <- list(control=list(dist_fn=dist, dist_fn_args=list(method="euclidean")),
#         lag=2,
#         theta=theta)
# 
#     expect_equal(ssr_predict(ssr_fit, train_data, predict_data, prediction_steps=prediction_step)$weights[[1]],
#         expected_weights)
# })
