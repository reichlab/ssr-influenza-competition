### numerical functions

## interface to R's C API for logspace arithmetic

logspace_sub <- function(logx, logy) {
	return(.Call("logspace_sub_C", as.numeric(logx), as.numeric(logy)))
}

logspace_add <- function(logx, logy) {
	return(.Call("logspace_add_C", as.numeric(logx), as.numeric(logy)))
}

logspace_sum_matrix_rows <- function(logX) {
	return(.Call("logspace_sum_matrix_rows_C", as.numeric(logX), as.integer(nrow(logX)), as.integer(ncol(logX))))
}

logspace_sub_matrix_rows <- function(logX) {
	if(!is.matrix(logX) || !identical(ncol(logX), 2L))
		stop("logX must be a matrix with 2 columns")

	return(.Call("logspace_sub_matrix_rows_C", as.numeric(logX), as.integer(nrow(logX))))
}


### forecast evaluation

#' Compute Mean Absolute Scaled Error (MASE) to evaluate point predictions
#' 
#' @param obs vector of observed values
#' @param pred vector of point predictions
mase <- function(obs, pred) {
    mean_forecast_error <- mean(abs(obs - pred))
    mean_naive_error <- mean(abs(obs[-length(obs)] - obs[-1]))
    return(mean_forecast_error / mean_naive_error)
}

#' Compute Mean Absolute Error to evaluate point predictions
#' 
#' @param obs vector of observed values
#' @param pred vector of point predictions
mae <- function(obs, pred) {
    if(length(obs) != length(pred)) {
        stop("obs and pred must have the same length")
    }
    
    return(mean(abs(obs - pred)))
}
