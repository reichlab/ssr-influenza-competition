### functions to fit ssr and make predictions and forecasts

### functions for creating parameters specifying how the fit should be performed

#' Assemble a list of ssr_control parameters for the ssr function with
#'     user-specified values.
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param max_lag a named list: names are the names of variables that can be
#'     used in the observation weighting process; entries give the maximum
#'     number of lag time steps that may included in the lagged observation
#'     vector for the corresponding variable
#' @param prediction_horizons integer vector: the number of time steps between
#'     the last observation and the time at which we make a prediction
#' @param kernel_fns a named list: names are the names of variables that can be
#'     used in the observation weighting process or the lead process; entries
#'     are the names of functions to use in computing the value of the kernel
#'     function for the corresponding variable (and its lagged versions).
#' @param theta_est a named list: one component for each variable in X_names,
#'     each component is a character vector with names of parameters to the
#'     corresponding kernel function that will be estimated
#' @param theta_fixed a named list: one component for each variable in X_names,
#'     each component is a named list with values of parameters for the
#'     corresponding kernel function that will be held fixed
#' @param theta_transform_fns a named list: one component for each kernel
#'     function that is used.  Each component is a named list of parameters for
#'     that kernel function which are being estimated.  Each component of that
#'     sub-list is a list with two components: "transform" is a function
#'     that transforms the given parameter to the estimation scale, and
#'     "detransform" is a function that transforms the parameter back to the
#'     original scale
#' @param crossval_buffer during cross-validation, the number of indices before
#'     the time at which we are making a prediction to drop from the "training
#'     examples".
#' @param loss_fn_name a string giving the name of the function use to
#'     compute loss from predictions
#' @param loss_fn_args a named list giving arguments to the loss function
#' @param par_packages a character vector containing names of packages that need
#'     to be loaded in instances of R when computations are performed in
#'     parallel.
#' 
#' @return the (at this point, unvalidated) list of ssr_control parameters
create_ssr_control <- function(X_names,
        y_names,
        time_name,
        max_lag,
        prediction_horizons,
        kernel_fns,
        theta_est,
        theta_fixed,
        theta_transform_fns,
        crossval_buffer,
        loss_fn_name,
        loss_fn_args,
        par_packages = NULL) {
    ssr_control <- list()

    ssr_control$X_names <- X_names
    ssr_control$y_names <- y_names
    ssr_control$time_name <- time_name
    
    ssr_control$max_lag <- max_lag
    ssr_control$prediction_horizons <- prediction_horizons
    
    ssr_control$kernel_fns <- kernel_fns
    ssr_control$theta_est <- theta_est
    ssr_control$theta_fixed <- theta_fixed
    ssr_control$theta_transform_fns <- theta_transform_fns
    
    ssr_control$crossval_buffer <- crossval_buffer
    
    ssr_control$loss_fn_name <- loss_fn_name
    ssr_control$loss_fn_args <- loss_fn_args

    return(ssr_control)
}

#' Assemble a list of ssr_control parameters for the ssr function with default
#'     values
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' 
#' @return the list of ssr_control parameters
create_ssr_control_default <- function(X_names, y_names, time_name, data) {
    ssr_control <- list()
    
    ssr_control$X_names <- X_names
    ssr_control$y_names <- y_names
    ssr_control$time_name <- time_name
    
    ssr_control$max_lag <- lapply(X_names, function(xn) 1)
    names(ssr_control$max_lag) <- X_names
    
    ssr_control$kernel_fns <- get_default_kernel_fns(X_names,
        y_names,
        time_name,
        data)
    
    ssr_control$theta_est <- lapply(ssr_control$kernel_fns, function(kernel_fn) {
        if(identical(kernel_fn, "squared_exp_kernel")) {
            return("bw")
        } else if(identical(kernel_fn, "discrete_kernel")) {
            return("bw")
        } else if(identical(kernel_fn, "periodic_kernel")) {
            return(c("bw", "period"))
        }
    })
    ssr_control$theta_fixed <- lapply(ssr_control$kernel_fns, function(kernel_fn) {
        return(NULL)
    })
    
    ssr_control$loss_fn_name <- "mase"
    ssr_control$loss_fn_args <- list()
    
    return(ssr_control)
}

#' Get default kernel functions based on a brief look at the data.  This is
#' unreliable.  Update to return periodic_kernel if X_names[i] == time_name?
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' 
#' @return the list of ssr_control parameters
get_default_kernel_fns <- function(X_names, y_names, time_name, data) {
    kernel_fns <- lapply(c(X_names, y_names), function(xyn) {
        if(isTRUE(all.equal(data[, xyn], as.integer(data[, xyn])))) {
            return("discrete_kernel")
        } else {
            return("squared_exp_kernel")
        }
    })
    names(kernel_fns) <- c(X_names, y_names)

    return(kernel_fns)
}

#' Validate ssr_control parameters for ssr -- not implemented
#' 
#' @param ssr_control a list of ssr_control parameters for ssr
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' 
#' @return no return value -- either stops with an error or not.
validate_ssr_control <- function(ssr_control, X_names, y_names, time_name, data) {
#    warning("ssr ssr_control parameter validation not yet implemented")
}


### functions for parameter estimation

#' Estimate the parameters for ssr.  There is redundancy here in that X_names,
#' y_names, time_name are all included in the ssr_control object as well as
#' parameters to this function.  Decide on an interface.
#' 
#' @param X_names a character vector of length >= 1 containing names of
#'     variables in the data data frame to use in forming the lagged
#'     observation process used for calculating weights
#' @param y_names a character vector of length 1 containing the name of the
#'     variable in the data data frame to use as the target for prediction
#' @param time_name (optional) a character vector of length 1 containing the
#'     name of the variable in the data data frame to use as the time.
#' @param data a data frame where rows are consecutive observations
#' @param ssr_control a list of parameters ssr_controlling how the fit is done.
#'     See the documentation for ssr_control.
#' 
#' @return an object representing an estimated ssr model.  Currently, a list
#'     with 7 components.
ssr <- function(X_names,
        y_names,
        time_name,
        data,
        ssr_control) {
    ## get/validate ssr_control argument
    if(missing(ssr_control)) {
        ssr_control <- ssr_control_default(X_names, y_names, time_name, data)
        warning("ssr_control argument not supplied to ssr -- using defaults, which may be bad")
    } else {
        validate_ssr_control(ssr_control, X_names, y_names, time_name, data)
    }
    
    ## estimate lags and kernel parameters via cross-validation
    param_estimates <- est_ssr_params_stepwise_crossval(data, ssr_control)
    
    return(list(ssr_control=ssr_control,
        X_names=X_names,
        y_names=y_names,
        time_name=time_name,
        lags_hat=param_estimates$lags_hat,
        theta_hat=param_estimates$theta_hat,
        train_data=data))
}

#' Use a stepwise procedure to estimate parameters and select lags by optimizing
#' a cross-validation estimate of predictive performance
#' 
#' @param data the data frame to use in performing cross validation
#' @param ssr_control a list of parameters specifying how the fitting is done
#' 
#' @return a list with two components: lags_hat is the estimated "optimal" lags
#'     to use for each variable, and theta_hat is the estimated "optimal"
#'     kernel parameters to use for each combination of variable and lag
est_ssr_params_stepwise_crossval <- function(data, ssr_control) {
    ## assemble estimation examples based on full data set with all possible
    ## lags -- this wastes memory if we have a lot of data, but saves
    ## us from having to re-calculate many times.
    all_lags_as_list <- lapply(ssr_control$X_names, function(xn) {
        lapply(seq(from=0, to=ssr_control$max_lag[[xn]]), function(lag_value) {
            list(var_name=xn,
                lag_value=lag_value)
        })
    }) %>% unlist(recursive=FALSE)
    
    ## initialize cross-validation process
    ## the variable selected in previous iteration
    selected_var_lag_ind <- NULL
    ## cross-validation estimate of loss associated with current estimates
    crossval_prediction_loss <- Inf
    ## initial parameters: no variables/lags selected, no kernel parameters
    lags_hat <- list()
    theta_hat <- list()
    
    repeat {
        ## get cross-validation estimates of performance for model obtained by
        ## adding or removing each variable/lag combination (except for the one
        ## updated in the previous iteration) from the model, as well as the
        ## corresponding parameter estimates
        
        ## commented out use of foreach for debugging purposes
        crossval_results <- foreach(i=seq_along(all_lags_as_list),
            .packages=c("ssr", ssr_control$par_packages),
            .combine="c") %dopar% {
#        crossval_results <- lapply(seq_along(all_lags_as_list), function(i) {
            if(length(selected_var_lag_ind) > 0 && i == selected_var_lag_ind) {
                list(
                    list(loss=crossval_prediction_loss,
                        lags=lags_hat,
                        theta=theta_hat)
                )
            } else {
               list( # put results in a list so that combine="c" is useful
                    est_ssr_params_stepwise_crossval_one_potential_step(
                        prev_lags=lags_hat,
                        prev_theta=theta_hat,
                        update_var_name=all_lags_as_list[[i]]$var_name,
                        update_lag_value=all_lags_as_list[[i]]$lag_value,
                        data=data,
                        ssr_control=ssr_control)
                )
            }
        }
#        })
        
        ## pull out loss achieved by each model, find the best value
        loss_achieved <- sapply(crossval_results, function(component) {
            component$loss
        })
        optimal_loss_ind <- which.min(loss_achieved)
        
#        print("loss achieved is:")
#        print(loss_achieved)
#        print("\n")
        
        ## either update the model and keep going or stop the search
        if(loss_achieved[optimal_loss_ind] < crossval_prediction_loss) {
            ## found a model improvement -- update and continue
            selected_var_lag_ind <- optimal_loss_ind
            crossval_prediction_loss <- loss_achieved[selected_var_lag_ind]
            lags_hat <- crossval_results[[selected_var_lag_ind]]$lags
            theta_hat <- crossval_results[[selected_var_lag_ind]]$theta
        } else {
            ## could not find a model improvement -- stop search
            break
        }
    }

    return(list(lags_hat=lags_hat,
        theta_hat=theta_hat))
}

#' Estimate the parameters theta and corresponding cross-validated estimate of
#' loss for one possible model obtained by adding or removing a variable/lag
#' combination from the model obtained in the previous iteration of the stepwise
#' search procedure.
#' 
#' @param prev_lags list representing combinations of variables and lags
#'     included in the model obtained at the previous step
#' @param prev_theta list representing the kernel parameter estimates obtained
#'     at the previous step
#' @param update_var_name the name of the variable to try adding or removing
#'     from the model
#' @param update_lag_value the value of the lag for the variable specified by
#'     update_var_name to try adding or removing from the model
#' @param data the data frame with observations used in estimating model
#'     parameters
#' @param ssr_control a list of parameters specifying how the fitting is done
#' 
#' @return a list with three components: loss is a cross-validation estimate of
#'     the loss associated with the estimated parameter values for the given
#'     model, lags is a list representing combinations of variables and lags
#'     included in the updated model, and theta is a list representing the
#'     kernel parameter estimates in the updated model
est_ssr_params_stepwise_crossval_one_potential_step <- function(prev_lags,
    prev_theta,
    update_var_name,
    update_lag_value,
    data,
    ssr_control) {
    
    updated_lags <- prev_lags
    updated_theta <- prev_theta
    
    update_var_lag_combo <- paste0(update_var_name, "_lag", update_lag_value)
    if(update_lag_value %in% prev_lags[[update_var_name]]) {
        ## remove variable/lag combination from model
        
        ## remove lag
        if(length(updated_lags[[update_var_name]]) == 1) {
            updated_lags <- updated_lags[
                -(names(updated_lags) == update_var_name)
            ]
        } else {
            updated_lags[[update_var_name]] <- updated_lags[[update_var_name]][
                updated_lags[[update_var_name]] != update_lag_value
            ]
        }
        
        ## remove parameter values
        updated_theta <- updated_theta[
            -(names(updated_theta) == update_var_lag_combo)
        ]
    } else {
        ## add variable/lag combination to model
        
        ## add lag
        updated_lags[[update_var_name]] <- sort(
            c(updated_lags[[update_var_name]], update_lag_value)
        )
        
        ## add default initial parameter values
        updated_theta[[update_var_lag_combo]] <-
            get_kernel_fn_init_params(update_var_name, ssr_control)
        
        ## ensure lags and theta are in order by variable name
        ## note that the lag values are not in the same order -- lags orders as
        ## numeric, but theta orders as character
        ## I don't actually rely on order right now, use names instead.
        ## Keeping this in case.
        updated_lags <- updated_lags[order(names(updated_lags))]
        updated_theta <- updated_theta[order(names(updated_theta))]
    }
    
    ## transform to estimation scale and vectorize for use by optim
    updated_theta <- transform_theta_list(updated_theta,
        updated_lags,
        ssr_control,
        direction="transform")
    theta_vector <- vectorize_theta(updated_theta)
    
    tryCatch({
        optim_result <- optim(par=theta_vector,
   	     fn=ssr_crossval_estimate_parameter_loss,
#	             gr = gradient_ssr_crossval_estimate_parameter_loss,
		gr=NULL,
		        lags=updated_lags,
        		data=data,
       			 ssr_control=ssr_control,
        		 method="L-BFGS-B",
        #		lower=-10000,
        #		upper=10000,
        #control=list(),
 	       hessian=FALSE)

        print(optim_result)
    
        ## convert back to list and original parameter scale
        updated_theta <- unvectorize_theta(optim_result$par,
            updated_lags,
            ssr_control,
            add_fixed_params=FALSE)
        updated_theta <- transform_theta_list(updated_theta,
            updated_lags,
            ssr_control,
            direction="detransform")
    
         return(list(
            loss=optim_result$value,
	        lags=updated_lags,
	        theta=updated_theta
	    ))
    }, error = function(e) {
       return(list(loss=Inf,
		lags=updated_lags,
		theta=updated_theta))
    })

    return(list(loss=Inf,
        lags=updated_lags,
	theta=updated_theta))
}

#' Using cross-validation, estimate the loss associated with a particular set
#' of lags and kernel function parameters.
#' 
#' @param theta_vector vector of kernel function parameters that are being
#'     estimated
#' @param lags list representing combinations of variables and lags
#'     included in the model
#' @param data the data frame to use in performing cross validation
#' @param ssr_control a list of parameters specifying how the fitting is done
#' 
#' @return numeric -- cross-validation estimate of loss associated with the
#'     specified parameters
ssr_crossval_estimate_parameter_loss <- function(theta_vector,
    lags,
    data,
    ssr_control) {
    ## set up theta list containing both the kernel parameters that are being
    ## estimated and the kernel parameters that are being held fixed
    ## also, transform back to original scale
    theta <- unvectorize_theta(theta_vector,
        lags,
        ssr_control,
        add_fixed_params=TRUE)
    theta <- transform_theta_list(theta,
        lags,
        ssr_control,
        direction="detransform")
    
    ## create data frame of "examples" -- lagged observation vectors and
    ## corresponding prediction targets
    cross_validation_examples <- assemble_training_examples(data,
        lags,
        ssr_control$y_names,
        leading_rows_to_drop=max(unlist(ssr_control$max_lag)),
        additional_rows_to_drop=NULL,
        prediction_horizons=ssr_control$prediction_horizons,
        drop_trailing_rows=TRUE)
    
    ## This could be made more computationally efficient by computing
    ## kernel values for all relevant combinations of lags for each variable,
    ## then combining as appropriate -- currently, the same kernel value may be
    ## computed multiple times in the call to ssr_predict_given_lagged_obs
    crossval_loss_by_time_ind <- sapply(
        seq_len(nrow(cross_validation_examples$lagged_obs)),
        function(t_pred) {
            ## get training indices -- those indices not within
            ## t_pred +/- ssr_control$crossval_buffer
            t_train <- seq_len(nrow(cross_validation_examples$lagged_obs))
            t_train <- t_train[!(t_train %in%
                seq(from=t_pred - ssr_control$crossval_buffer,
                    to=t_pred + ssr_control$crossval_buffer))]
            
            ## calculate kernel weights and centers for prediction at
            ## prediction_lagged_obs based on train_lagged_obs and
            ## train_lead_obs
            ## assemble lagged and lead observations -- subsets of
            ## cross_validation_examples given by t_pred and t_train
            ## we can re-use the weights at different prediction_target_inds,
            ## and just have to adjust the kernel centers
            prediction_target_ind <- 1
            
            train_lagged_obs <- cross_validation_examples$lagged_obs[
                t_train, , drop=FALSE]
            train_lead_obs <- cross_validation_examples$lead_obs[
                t_train, prediction_target_ind, drop=FALSE]
            prediction_lagged_obs <- 
                cross_validation_examples$lagged_obs[
                    t_pred, , drop=FALSE]
            
            kernel_weights_and_centers <- ssr_predict_given_lagged_obs(
                train_lagged_obs=train_lagged_obs,
                train_lead_obs=train_lead_obs,
                prediction_lagged_obs=prediction_lagged_obs,
                ssr_fit=list(lags_hat=lags,
                    theta_hat=theta,
                    ssr_control=ssr_control
                ),
                normalize_weights=TRUE)
            
            ## for each prediction target variable, compute loss
            crossval_loss_by_prediction_target <- sapply(
                seq_len(ncol(cross_validation_examples$lead_obs)),
                function(prediction_target_ind) {
                    ## update kernel centers to match prediction_target_ind
                    kernel_weights_and_centers$centers <-
                        cross_validation_examples$lead_obs[
                            t_train, prediction_target_ind, drop=FALSE]
                    
                    ## calculate and return value of loss function based on prediction
                    ## and realized value
                    loss_fn_args <- ssr_control$loss_fn_args
                    loss_fn_args$kernel_weights_and_centers <- list(
                        weights=kernel_weights_and_centers$weights,
                        centers=kernel_weights_and_centers$centers[, 1]
                    )
                    
                    loss_fn_args$obs <- as.numeric(
                        cross_validation_examples$lead_obs[
                            t_pred, prediction_target_ind]
                    )
                    
                    return(do.call(ssr_control$loss_fn_name, loss_fn_args))
                })
            
            return(sum(crossval_loss_by_prediction_target))
        })
    
#    browser()
    cat(sum(crossval_loss_by_time_ind))
    cat(is.na(sum(crossval_loss_by_time_ind)))
    cat("\n")
    
    if(is.finite(sum(crossval_loss_by_time_ind)) && !is.na(sum(crossval_loss_by_time_ind))) {
        return(sum(crossval_loss_by_time_ind))
    } else {
      	cat("returning double max: ")
	cat(.Machine$double.xmax)
	cat("\n")
        return(.Machine$double.xmax)
    }
}

#' Get initial parameter values in list form for the given kernel function
#' 
#' @param var_name the name of the variable for which we are getting kernel
#'     function parameters
#' @param ssr_control control parameters for the ssr fit
#' 
#' @param list of parameter values for the given kernel function
get_kernel_fn_init_params <- function(var_name, ssr_control) {
    return_val <- list()
    for(i in seq_along(ssr_control$theta_est[[var_name]])) {
        return_val[[ssr_control$theta_est[[var_name]][i]]] <- 0.1
    }
    
    return(return_val)
}

#' Convert theta from list form to vector form.
#' 
#' @param theta_list kernel parameters theta in list form
#' 
#' @return numeric vector with parameter values
vectorize_theta <- function(theta_list) {
    return(unlist(theta_list))
}

#' Convert theta from vector form to list form
#' 
#' @param theta_vector kernel parameters theta in vector form
#' @param lags list representing combinations of variables and lags
#'     included in the model
#' @param ssr_control control parameters for the ssr fit
#' @param add_fixed_params boolean -- should parameters that are being held
#'     fixed in the estimation process be added to the return value?
#' 
#' @return list of lists of parameter values -- outer list has one component
#'     for each combination of variable and lag, inner list has one component
#'     for each parameter used in the corresponding kernel function
unvectorize_theta <- function(theta_vector,
    lags,
    ssr_control,
    add_fixed_params) {
    theta_names <- rbind.fill(lapply(seq_along(lags), function(ind) {
        data.frame(var_name=names(lags)[ind],
            lag_val=lags[[ind]],
            combined_name=paste0(names(lags)[ind], "_lag", lags[[ind]]))
    }))
    
    theta_list <- lapply(seq_len(nrow(theta_names)), function(i) {})
    
    for(theta_list_ind in seq_len(nrow(theta_names))) {
        ## parameters that are being estimated
        theta_list[[theta_list_ind]] <-
            unvectorize_theta_one_kernel_fn(theta_vector,
                theta_names[theta_list_ind, "var_name"],
                theta_names[theta_list_ind, "lag_val"],
                ssr_control)

        ## parameters that are being held fixed, if requested
        if(add_fixed_params) {
            theta_list[[theta_list_ind]] <- c(theta_list[[theta_list_ind]],
                ssr_control$theta_fixed[[
                    as.character(theta_names[theta_list_ind, "var_name"])
                ]]
            )
        }
    }
    
    names(theta_list) <- theta_names$combined_name

    return(theta_list)
}

transform_theta_list <- function(theta_list,
    lags,
    ssr_control,
    direction="transform") {
    ## all names of elements in theta_list, split into variable name and lag val
    ## easier and more reliable to do this based on lags
    theta_names <- rbind.fill(lapply(seq_along(lags), function(ind) {
        data.frame(var_name=names(lags)[ind],
            lag_val=lags[[ind]],
            combined_name=paste0(names(lags)[ind], "_lag", lags[[ind]]))
    }))
    
    ## obtain transformed list of parameters
    transformed_theta_list <- lapply(seq_len(nrow(theta_names)),
        function(component_ind) {
            ## original list of parameter values
            component_params <- theta_list[[component_ind]]
            ## name of variable used
            component_var <- theta_names$var_name[component_ind]
            ## name of kernel function used
            kernel_fn_name <- ssr_control$kernel_fns[[component_var]]
            
            transformed_component_params <- component_params
            for(param_name in names(component_params)) {
                ## transformation only necessary for parameters being estimated
                if(param_name %in% ssr_control$theta_est[[component_var]]) {
                    transform_fn_name <- ssr_control$theta_transform_fns[[
                        kernel_fn_name]][[param_name]][[direction]]
                    transformed_component_params[[param_name]] <-
                        do.call(transform_fn_name,
                            list(x=component_params[[param_name]]))
                }
            }
            return(transformed_component_params)
        }
    )
    
    names(transformed_theta_list) <- theta_names$combined_name
    
    return(transformed_theta_list)
}

#' Convert theta from vector form to list form for one kernel function
#' 
#' @param theta_vector kernel parameters theta in vector form for all kernel
#'     functions
#' @param theta_vector_ind the index in theta_vector to start getting parameters
#'     for the given kernel function
#' @param var_name the name of the variable for which we are getting kernel
#'     function parameters
#' @param lag_value the lag value for which we are getting kernel function
#'     parameters
#' @param ssr_control control parameters for the ssr fit
#' 
#' @return list of parameter values, one named component
#'     for each parameter used in the given kernel function
unvectorize_theta_one_kernel_fn <- function(theta_vector,
    var_name,
    lag_value,
    ssr_control) {
    return_val <- list()
    
    ## iterate over the parameters that are being estimated
    for(i in seq_along(ssr_control$theta_est[[var_name]])) {
        ## name of parameter
        param_name <- ssr_control$theta_est[[var_name]][i]
        ## combined variable name, lag value, and parameter name
        ## used as lookup name in vector
        vector_elt_name <- paste0(var_name, "_lag", lag_value, ".", param_name)
        ## store parameter value
        return_val[[param_name]] <- unname(theta_vector[vector_elt_name])
    }
    
    return(return_val)
}


#' Make predictions from an estimated ssr model forward prediction_horizon time
#' steps from the end of predict_data, based on the weighting variables, lags,
#' kernel functions, and bandwidths specified in the ssr_fit object.
#' 
#' @param ssr_fit is an object representing a fitted ssr model
#' @param prediction_data is a vector of data points to use in prediction
#' @param prediction_horizon is an integer specifying the number of steps ahead
#'     to perform prediction
#' @param normalize_weights boolean, should the weights be normalized?
#' 
#' @return a list with three components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
ssr_predict <- function(ssr_fit,
        prediction_data,
        leading_rows_to_drop=max(unlist(ssr_fit$lags_hat)),
        additional_training_rows_to_drop=NULL,
        prediction_horizon,
        normalize_weights=TRUE) {
    ## get training and prediction examples
    training_examples <- assemble_training_examples(ssr_fit$train_data,
        ssr_fit$lags_hat,
        ssr_fit$y_names,
        leading_rows_to_drop,
        additional_training_rows_to_drop,
        prediction_horizon,
        drop_trailing_rows=TRUE)
    
    prediction_examples <- assemble_prediction_examples(prediction_data,
        ssr_fit$lags_hat,
        leading_rows_to_drop)
    
    ## assemble kernel centers and weights
    ssr_predict_given_lagged_obs(training_examples$lagged_obs,
        training_examples$lead_obs,
        prediction_examples$lagged_obs,
        ssr_fit,
        prediction_horizon)
}

#' Construct data frame of lagged observation vectors and data frame of
#' corresponding prediction targets.
#' 
#' @param data data frame with observations of all variables at all times
#' @param lags list representing combinations of variables and lags
#'     included in the lagged observation vectors
#' @param y_names names of variables included as prediction targets
#' @param leading_rows_to_drop an integer specifying the number of leading rows
#'     to drop.  This will typically be max(unlist(lags)), but could be larger
#'     for example if we are in the process of searching lags to determine
#'     which ones to include.  Then for example if our search is for lags up to
#'     52, we would want to set leading_rows_to_drop = 52.
#' @param additional_rows_to_drop an integer vector specifying indices of
#'     additional rows to drop.  For example, if we are performing
#'     cross-validation, we might want to drop all rows within +/- 52 indices of
#'     the current prediction target.
#' @param prediction_horizons an integer vector specifying the number of time
#'     steps between the last observation and the prediction target.
#' @param drop_trailing_rows boolean:  drop the last prediction_horizon rows?
#'     These are the rows for which we can form lagged observation vectors, but
#'     we cannot obtain a corresponding prediction target.
#'     
#' @return a list with two components: lagged_obs is a data frame with lagged
#'     observation vectors, and lead_obs is a data frame with corresponding
#'     prediction targets.
assemble_training_examples <- function(data,
    lags,
    y_names,
    leading_rows_to_drop,
    additional_rows_to_drop,
    prediction_horizons,
    drop_trailing_rows=TRUE) {
    ## which rows should not be used as regression/density estimation examples
    ## either because the corresponding regression example cannot be formed
    ## or because we are performing cross-validation and don't want to use
    ## times adjacent to the prediction target

    ## too early
    all_train_rows_to_drop <- seq_len(leading_rows_to_drop)

    ## passed in indices -- near prediction target
    all_train_rows_to_drop <- c(all_train_rows_to_drop, additional_rows_to_drop)
    
    ## too late
    max_prediction_horizon <- max(prediction_horizons)
    if(drop_trailing_rows) {
        all_train_rows_to_drop <- c(all_train_rows_to_drop,
            seq(from=nrow(data) - max_prediction_horizon + 1,
                to=nrow(data)))
    }
    
    ## compute lagged observation vectors for train and prediction data
    train_lagged_obs <- compute_lagged_obs_vecs(data,
        lags,
        all_train_rows_to_drop)
    
    ## compute lead observation series for train data
    train_lead_obs <- data.frame(rep(NA, nrow(data)))
    for(y_name in y_names) {
        for(prediction_horizon in prediction_horizons) {
            train_lead_obs_inds <-
                seq(from=1, to=nrow(data)) + prediction_horizon
            component_name <- paste0(y_name, "_horizon", prediction_horizon)
            train_lead_obs[[component_name]] <-
                data[train_lead_obs_inds, y_name, drop=TRUE]
        }
    }
    train_lead_obs <- train_lead_obs[-all_train_rows_to_drop, , drop=FALSE]
    
    ## drop initial column of NA's
    train_lead_obs <- train_lead_obs[, -1, drop=FALSE]
    
    return(list(lagged_obs=train_lagged_obs,
        lead_obs=train_lead_obs))
}

#' Construct data frame of lagged observation vectors.
#' 
#' @param data data frame with observations of all variables at all times
#' @param lags list representing combinations of variables and lags
#'     included in the lagged observation vectors
#' @param leading_rows_to_drop an integer specifying the number of leading rows
#'     to drop.  This will typically be max(unlist(lags)), but could be larger
#'     for example if we are in the process of searching lags to determine
#'     which ones to include.  Then for example if our search is for lags up to
#'     52, we would want to set leading_rows_to_drop = 52.
#'     
#' @return a list with one component: lagged_obs is a data frame with lagged
#'     observation vectors.
assemble_prediction_examples <- function(data,
    lags,
    leading_rows_to_drop) {
    all_prediction_rows_to_drop <- seq_len(leading_rows_to_drop) # too early

    prediction_lagged_obs <- compute_lagged_obs_vecs(data,
        lags,
        all_prediction_rows_to_drop)
    
    return(list(lagged_obs=prediction_lagged_obs))
}

#' Make predictions from an estimated ssr model forward prediction_horizon time
#' steps from the end of predict_data, based on the kernel functions and
#' bandwidths specified in the ssr_fit object.  This function requires that the
#' lagged and lead observation vectors have already been computed.
#' 
#' @param train_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the training data.  Each row
#'     corresponds to a time point.  Each column is a (lagged) variable.
#' @param train_lead_obs is a vector with length = nrow(train_lagged_obs) with
#'     the value of the prediction target variable corresponding to each row in
#'     the train_lagged_obs matrix.
#' @param prediction_lagged_obs is a matrix (with column names) containing the
#'     lagged observation vector computed from the prediction data.  There is
#'     only one row, representing one time point.  Each column is a (lagged)
#'     variable.
#' @param ssr_fit is an object representing a fitted ssr model
#' @param normalize_weights boolean, should the weights be normalized?
#' 
#' @return a list with three components:
#'     log_weights: a vector of length = length(train_lagged_obs) with
#'         the log of weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     weights: a vector of length = length(train_lagged_obs) with
#'         the weights assigned to each observation (up to a constant of
#'         proportionality if normalize_weights is FALSE)
#'     centers: a copy of the train_lead_obs argument -- kernel centers
ssr_predict_given_lagged_obs <- function(train_lagged_obs,
    train_lead_obs,
    prediction_lagged_obs,
    ssr_fit,
    normalize_weights=TRUE) {
    
    ## compute log of kernel function values representing similarity
    ## of lagged observation vectors from train_data and predict_data.
    ## result is a vector of length nrow(train_lagged_obs)
    log_weights <- compute_kernel_values(train_lagged_obs,
        prediction_lagged_obs,
        ssr_fit$lags_hat,
        ssr_fit$theta_hat,
        ssr_fit$ssr_control,
        log = TRUE)

    ## if requested, normalize log weights
    if(normalize_weights) {
        log_weights <- compute_normalized_log_weights(log_weights)
    }
    
    return(list(log_weights=log_weights,
        weights=exp(log_weights),
        centers=train_lead_obs))
}

#' Normalize a vector of weights on log scale: given a vector of log_weights
#' such that exp(log_weights) are proportional to the final weights, update
#' so that exp(log_weights) sums to 1.
#' 
#' @param log_weights: a vector of log(w) where w is proportional to the weights
#' 
#' @return normalized log_weights so that sum(exp(log_weights)) = 1
compute_normalized_log_weights <- function(log_weights) {
    norm_const <- logspace_sum_matrix_rows(matrix(log_weights, nrow = 1))
    return(log_weights - norm_const)
}

#' Compute kernel values measuring the similarity of each row in the
#' train_lagged_obs data frame to the prediction_lagged_obs.
#' 
#' @param train_lagged_obs a data frame with lagged observation vectors computed
#'     from the training data
#' @param prediction_lagged_obs a data frame with the lagged observation vector
#'     computed from the prediction data.  It is assumed that
#'     prediction_lagged_obs contains only one row.
#' @param theta a named list with one component for each entry of lags
#'     (combined variable name and lag size).  This component is a named list
#'     with arguments to the kernel function.
#' @param ssr_control a list of ssr_control parameters for ssr
#' @param log boolean; if TRUE (default), return kernel values on the log scale
compute_kernel_values <- function(train_lagged_obs,
    prediction_lagged_obs,
    lags,
    theta,
    ssr_control,
    log = TRUE) {
    ## create a matrix of log kernel values by component
    log_kernel_component_values <- matrix(NA,
        nrow=nrow(train_lagged_obs),
        ncol=ncol(train_lagged_obs))
    colnames(log_kernel_component_values) <- names(train_lagged_obs)
    
    l_ply(seq_along(lags), function(ind) {
        vname <- names(lags)[ind]
        l_ply(lags[[ind]], function(lag_val) {
            col_name <- paste0(vname, "_lag", lag_val)
            
            ## assemble arguments to kernel function
            kernel_fn_args <- theta[[col_name]]
            kernel_fn_args$x <- train_lagged_obs[[col_name]]
            kernel_fn_args$center <- prediction_lagged_obs[[col_name]]
            kernel_fn_args$log <- TRUE
            
            ## call kernel function and store results
            log_kernel_component_values[, col_name] <<- 
                do.call(ssr_control$kernel_fns[[vname]], kernel_fn_args)
        })
    })
    
    ## return on scale requested by user
    ## these computations assume product kernel --
    ## if we're doing something else, change apply(..., sum)
    if(log) {
        return(apply(log_kernel_component_values, 1, sum))
    } else {
        return(exp(apply(log_kernel_component_values, 1, sum)))
    }
}

#' Compute a data frame with lagged observation vectors
#' 
#' @param data a data frame
#' @param lags a named list: The component name matches the name of 
#'     one of the variables in data, and the component value is an integer
#'     vector of lags to include for that variable
#' @param rows_to_drop an integer vector specifying rows to drop after computing
#'     lagged observation vectors.
compute_lagged_obs_vecs <- function(data,
    lags,
    rows_to_drop) {
    ## set or validate leading_rows_to_drop
    unlisted_lags <- unlist(lags)
    max_entry_in_lags <- max(unlisted_lags)
    if(any(rows_to_drop < 0 | rows_to_drop > nrow(data))) {
        stop("all entries of rows_to_drop must integers between 1 and nrow(data)")
    } else if(!all(seq_len(max_entry_in_lags) %in% rows_to_drop)) {
        stop("all integers between 1 and the maximum entry in the lags argument must be contained in rows_to_drop")
    }
    
    ## create a data frame with one column for each entry in the lags argument
    result_names <- lapply(seq_along(lags), function(ind) {
        paste0(names(lags)[ind], "_lag", lags[[ind]])
    }) %>% unlist
    result <- as.data.frame(matrix(NA,
        nrow=nrow(data),
        ncol=length(result_names)))
    names(result) <- result_names
    
    ## set column values in result
    l_ply(seq_along(lags), function(ind) {
        vname <- names(lags)[ind]
        l_ply(lags[[ind]], function(lag_val) {
            result_name <- paste0(vname, "_lag", lag_val)
            result[seq(from=lag_val + 1, to=nrow(result)), result_name] <<-
                data[seq(from=1, to=nrow(result) - lag_val), vname]
        })
    })
    
    ## drop specified rows
    result <- result[-rows_to_drop, , drop=FALSE]
    
    return(result)
}

#' Squared exponential kernel function
#' 
#' @param x a vector of values at which to evaluate the kernel function
#' @param center a real number, center point for the kernel function
#' @param bw kernel bandwidth
#' @param return log scale value?
#' 
#' @return vector of the kernel function value at each point in x
squared_exp_kernel <- function(x, center, bw, log) {
    if(!is.numeric(center) | length(center) != 1) {
        stop("center must be numeric with length 1")
    }
    
    result <- -0.5 * ((x - center) / bw)^2
    
    if(log) {
        return(result)
    } else {
        return(exp(result))
    }
}

discrete_kernel <- squared_exp_kernel # replace with something real later

#' Periodic kernel function
#' 
#' @param x a vector of values at which to evaluate the kernel function
#' @param center a real number, center point for the kernel function
#' @param period kernel period
#' @param bw kernel bandwidth
#' 
#' @return vector of the kernel function value at each point in x
periodic_kernel <- function(x, center, period, bw, log) {
    result <- -0.5 * (sin(period * (x - center)) / bw)^2
    
    if(log) {
        return(result)
    } else {
        return(exp(result))
    }
}

#' This is a wrapper for the ssr_predict function to perform prediction from
#' the dengue data sets for the competition.  Computes KDE predictions for
#' the number of cases in the weeks indexed by t + prediction_horizon, where
#'   - t is specified by last_obs_season and last_obs_week
#'   - prediction_horizon varies over the values in prediction_horizons
#' Currently the function assumes that data have already been smoothed on the
#' log scale and a data frame variable named smooth_log_cases is available.
#' Log cases are used to do the SSR and get weights for each time point.
#' KDE estimates are then computed both on the log scale using smooth_log_cases
#' and also on the original data scale using exp(smooth_log_cases) as the
#' centers.
#' 
#' @param last_obs_season is the season for the last observed week before we
#'     want to make a prediction.  Has the form "2008/2009"
#' @param last_obs_week is the last observed week of the season specified by
#'     last_obs_season before we want to make a prediction.  An integer
#'     between 1 and 52.
#' @param lag_max is the number of lags to use in forming the state space
#'     representation via lagged observations
#' @param prediction_horizons is a vector giving the number of steps ahead
#'     to do prediction from the week before first_predict_week
#' @param data is the data set to use
#' @param season_var is a string with the name of the variable in the data data
#'     frame with the season.  I'm sure there is a better way to do this...
#' @param week_var is a string with the name of the variable in the data data
#'     frame with the week in the season.  I'm sure there is a better way to do this...
#' @param predict_vars is a string vector with the names of the variables in the data data
#'     frame to use for SSR.  I'm sure there is a better way to do this...
#' @param prediction_types is a character vector indicating the prediction types
#'     to perform: may contain one or more of "pt" and "density"
#' 
#' @return list of data frames with predictions.  one list component per entry
#'     in prediction_types argument.  Density estimates are from KDE along a
#'     grid of values for total_counts for each week that we predicted at.
ssr_predict_dengue_one_week <- function(last_obs_season,
    last_obs_week,
    lags,
    theta,
    prediction_horizon,
    data,
    season_var,
    week_var,
    X_names,
    y_names,
    time_name,
    prediction_types = c("pt", "density")) {
    
  	if(length(y_names) > 1) {
  		stop("SSR Prediction for multiple variables is not yet implemented!")
  	}
	
    ## get indices for prediction data
    ## in order to predict at time t + prediction_horizon,
    ##  - data must start at time t - lag_max so that we can get lag_max + 1
    ##    values in forming the lagged observation vector
    ##  - we need the number of data points equal to lag_max + 1
    ## first_predict_ind is the first index in the data data frame
    ## that needs to be included in the prediction data
    lag_max <- max(unlist(lags))
    last_obs_ind <- which(data[[season_var]] == last_obs_season &
        data[[week_var]] == last_obs_week)
    first_predict_data_ind <- last_obs_ind - lag_max
    predict_data_inds <- seq(from=first_predict_data_ind, length=lag_max + 1)
    
    ## for training inds, drop anything within +/- 1 year of the last
    ## observation
    time_points_per_year <- 52
    predict_seasons <- data[[season_var]][last_obs_ind + prediction_horizon]
    additional_training_rows_to_drop <-
        seq(from=max(1, last_obs_ind - time_points_per_year),
            to=min(nrow(data), last_obs_ind + time_points_per_year))
    
    ## do prediction
    ssr_predictions <- ssr_predict(
        ssr_fit=list(ssr_control=ssr_control_default(X_names, y_names, time_name, data),
            X_names=X_names,
            y_names=y_names,
            time_name=time_name,
            lags_hat=lags,
            theta_hat=theta,
            train_data=data),
        prediction_data=data[predict_data_inds, ],
        leading_rows_to_drop=max(unlist(lags)),
        additional_training_rows_to_drop=additional_training_rows_to_drop,
        prediction_horizon=prediction_horizon,
        normalize_weights=TRUE)
    
    ## get point estimate or density predictions
    result <- list()
    
    if("pt" %in% prediction_types) {
        result$pt_preds <- get_pt_predictions_one_week(ssr_predictions)
    }
    
    if("density" %in% prediction_types) {
        result$dist_preds <- get_dist_predictions_one_week(ssr_predictions)
    }
    
    return(result)
}

## a function to get point predictions for one week
get_pt_predictions_one_week <- function(ssr_predictions) {
    # get weighted mean
    pt_est <- weighted.mean(ssr_predictions$centers,
        ssr_predictions$weights)

    return(pt_est)
}

## a function to get kde predictions for one week
get_dist_predictions_one_week <- function(ssr_predictions) {
    ## do weighted kde
    kde_est <- density(ssr_predictions$centers,
        weights = ssr_predictions$weights,
        bw = "SJ")

    return(data.frame(x = kde_est$x,
        est_density = kde_est$y,
        est_bw=kde_est$bw))
}


#' Get week and season when we're predicting based on the index of the last obs
#' and the number of steps forward we're predicting
#' 
#' @param last_obs_ind index in the data data frame of the last observation
#'     before we start predicting
#' @param prediction_horizon the number of steps forward we're predicting
#' @param data the data frame, with columns named season_week and season
#' @param season_var is a string with the name of the variable in the data data
#'     frame with the season.  I'm sure there is a better way to do this...
#' @param week_var is a string with the name of the variable in the data data
#'     frame with the week in the season.  I'm sure there is a better way to do this...
#' 
#' @return a list with two components indicating the time of prediction:
#'     1) week is an integer from 1 to 52 and
#'     2) season is a string of the form "2008/2009"
get_prediction_season_week <- function(last_obs_ind, prediction_horizon, data, season_var, week_var) {
    wk <- (data[[week_var]][last_obs_ind] + prediction_horizon) %% 52
    if(wk == 0) {
        wk <- 52
    }
    seasons_advanced <- (data[[week_var]][last_obs_ind] + prediction_horizon
            - wk) / 52
    start_season_last_obs <- as.integer(substr(
            as.character(data[[season_var]][last_obs_ind]),
            start = 1,
            stop = 4))
	season <- start_season_last_obs + seasons_advanced
	
	if(nchar(as.character(data[[season_var]][last_obs_ind])) > 4) {
		season <- paste0(season, "/", season + 1)
	}
    
    return(list(week = wk, season = season))
}


### loss functions

#' Compute MASE from predictions that are in the form of kernel weights and
#' centers.  This is currently broken because we need predictions for a whole
#' time series to compute MASE, and the current "implementation" only takes
#' values for one time point....
mase_from_kernel_weights_and_centers <- function(
    kernel_weights_and_centers,
    obs) {
    stop("mase_from_kernel_weights_and_centers is broken -- do not use")
    pred <- get_pt_predictions_one_week(kernel_weights_and_centers)
    mase(obs, pred)
}

#' Compute "mean" absolute error for one time point from prediction in the form
#' of kernel weights and centers.
#' 
#' @param kernel_weights_and_centers a named list with two components: weights
#'     is a vector of kernel weights, and centers is a vector of kernel centers
#' @param obs is a numeric with length 1 containign the observed value for one
#'     time point
#' 
#' @return abs(obs - prediction) where prediction is the weighted mean of the
#'     kernel centers.
mae_from_kernel_weights_and_centers <- function(
    kernel_weights_and_centers,
    obs) {
    pred <- get_pt_predictions_one_week(kernel_weights_and_centers)
    mae(obs, pred)
}

#' Get the indices of the smallest k elements of v.  This code currently assumes
#' that k >= length(v)
#' 
#' @param v a vector
#' @param k number of indices to return
#' 
#' @return a vector of length k containing the indices of the k smallest
#'   elements of v, in ascending order.
get_inds_smallest_k <- function(v, k) {
    return(order(v, decreasing=FALSE)[seq_len(k)])
}
