#' Run all analyses
#' 
#' By default will run analyses on all datasets. It expects an input data frame with columns of time and rate, each row representing a single observation. To handle multiple datasets at once, have a citation column indicating the name for the dataset (which could be a full citation or as simple as "dataset1", "dataset2", etc.). You can also include a column for numerator (the number of events) and denominator (the number of opportunities for events). If you have a column for total_time, it will use that for randomization rather than time (this is important for things like phylogenetic trees, where the time over which events happen is not just the overall age of the tree but the sum of its branch lengths). 
#' 
#' Data can be randomized to help suggest whether an observed pattern is spurious. 
#' @param all_data A data frame with columns of time, rate, and perhaps citation, numerator, denominator, and/or total_time (other columns will be ignored).
#' @param nreps The number of times to randomize the data within each dataset.
#' @param epsilon_lower The lower bound for the epsilon parameter.
#' @return A data.frame of results with class hyperr8.
#' @export
#' @examples
#' library(hyperr8)
#' car_data <- generate_car_simulation()
#' all_run <- hyperr8_run(car_data)
#' plot(all_run)
hyperr8_run <- function(all_data, nreps=5, epsilon_lower=0) {
	all_data <- clean_input_data(all_data)
	original <- summarize_all_fitted_models(optimization_over_all_data(all_data, epsilon_lower=epsilon_lower))
	original$rep <- "Original"
	randomized <- optimization_and_summarization_over_randomized_data(all_data, nreps=nreps, epsilon_lower=epsilon_lower)
	randomized$rep <- paste0("Rep ", randomized$rep)
	merged <- dplyr::bind_rows(original, randomized)
	class(merged) <- c("hyperr8", class(merged))
	#merged <- as.data.frame(merged) # I am a Klingon with bad spelling. Death to tibbles. They have no honor.
	return(merged)
}

get_best_empirical <- function(run_output) {
	return(subset(run_output, rep=="Original" & deltaAIC==0))
}

#' Plotting function for hyperr8
#' 
#' This will plot the results from hyperr8_run.
#' @param x The output from hyperr8_run.
#' @param loglog Whether to use a log-log plot
#' @param ... Other arguments to pass to the plotting function.
#' @return A ggplot2 object.
#' @export
plot.hyperr8 <- function(x, loglog=TRUE,...) {
	npoints <- nrow(subset(x, rate_type=='empirical_rate' & deltaAIC==0))
	alpha <- min(0.2, 10/sqrt(npoints))
	gcool <- ggplot(subset(x, rate_type=='empirical_rate' & deltaAIC==0), aes(x=time, y=rate)) + geom_point(alpha=alpha) + facet_grid(dataset~rep) + theme_bw() + xlab("Time") + ylab("Rate")
	if(loglog) {
		gcool <- gcool + scale_x_continuous(trans = "log") + scale_y_continuous(trans = "log")
	}
	gcool <- gcool + geom_line(data=subset(x, rate_type=='predicted_rate' & deltaAIC==0), aes(x=time, y=rate, group=rep, colour=model))
	return(gcool)
}

#' Car simulation
#'
#' This will generate a data frame with simulated car data of speed, distance, and time. 
#' By default, it will have 1000 rows representing 1000 cars. Each car will have a simulated driving time and a simulated driving speed. The driving time and driving speed will be correlated. The correlation will be 0.99 by default.
#' @param mean_driving_time The mean driving time for all cars.
#' @param mean_driving_speed The mean driving speed for all cars.
#' @param sd_driving_time The standard deviation of the driving time for the cars.
#' @param sd_driving_speed The standard deviation of the driving speed for the cars.
#' @param cor_driving The correlation between driving time and driving speed.
#' @param n The number of cars to simulate.
#' @return A data frame with simulated car data.
#' 
#' @examples
#' library(hyperr8)
#' car_data <- generate_car_simulation()
#' 
#' # Original data
#' plot(car_data$time, car_data$distance, pch=".", col=rgb(0,0,0,0.2))
#' 
#' # Estimated speed
#' plot(car_data$time, car_data$rate, pch=".", col=rgb(0,0,0,0.2))
#' 
#' # Log-log estimated speed
#' plot(car_data$time, car_data$rate, pch=".", col=rgb(0,0,0,0.2), log="xy")
#' @export 
generate_car_simulation <- function(mean_driving_time=3, mean_driving_speed=70, sd_driving_time=0.2, sd_driving_speed=30, cor_driving=0.99, n=1000) {
	mean_driving_distance <- mean_driving_speed*mean_driving_time
	sd_driving_time <- 0.2
	sd_driving_distance <- 30
	cor_driving <- 0.99

	Sigma <- matrix(c(sd_driving_distance, cor_driving, cor_driving, sd_driving_time), nrow=2)
	mu <- c(mean_driving_distance, mean_driving_time)
	n <- 1000
	sims <- as.data.frame(MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma, empirical=TRUE))
	generating_speed<- mean_driving_distance/mean_driving_time
	names(sims) <- c("distance", "time")
	sims$rate <- sims$distance/sims$time
	sims$citation <- "simulated_car"
	return(sims)
}

# #former f1
# function_constant <- function(par, focal_data, do_log1p=TRUE) {
# 	varepsilon_0_plus_k <- par[1]
# 	#log(\hat{r}(t)) = log(\varepsilon_0 + k) - log(\hat{t})
# 	if(do_log1p) {
# 		return(log(varepsilon_0_plus_k) - focal_data$log_time)
# 	} else {
# 		return(varepsilon_0_plus_k/exp(focal_data$log_time))
# 	}
# }


	
# #former f2
# function_hyperbola <- function(par, focal_data, do_log1p=TRUE) {
# 	varepsilon_0 <- par[1]
# 	k <- par[2]
# 	if(do_log1p) {
# 		return(log(varepsilon_0 + k/exp(focal_data$log_time)) - focal_data$log_time)
# 	} else {
# 		return((varepsilon_0 + k/exp(focal_data$log_time))/exp(focal_data$log_time))
# 	}
# }

# #former f3
# function_linear <- function(par, focal_data, do_log1p=TRUE) {
# 	varepsilon_0 <- par[1]
# 	k <- par[2]
# 	if(do_log1p) {
# 		return(log(varepsilon_0 + k*exp(focal_data$log_time)) - focal_data$log_time)
# 	} else {
# 		return((varepsilon_0 + k*exp(focal_data$log_time))/exp(focal_data$log_time))
# 	}
# }


#' Function to find a constant for which log(y+c) is as symmetric as possible at the alpha and 1-alpha percentiles
#' 
#' This is done rather than a fixed log1p offset or a fixed tiny offset.
#' Code is from:
#' whuber (https://stats.stackexchange.com/users/919/whuber), Choosing c such that log(x + c) would remove skew from the population, URL (version: 2012-11-05): https://stats.stackexchange.com/q/41415
#' @param y A vector of values
#' @param alpha The percentile to use for the lower bound
#' @param beta The percentile to use for extrapoloation
#' @param gamma The amount to extrapolate below the min
#' @return A list with the offset and a code

log.start <- function(y, alpha=0.25, beta=1, gamma=0.01) {
  #
  # Find a constant `const` for which log(y+c) is as symmetric
  # as possible at the `alpha` and 1-alpha percentiles.
  # Returns a list containing the solution and a code:
  #   Code=0: Constant is good.
  #   Code=1: Data are already symmetric.
  #   Code=2: The best solution is not large enough and has been adjusted.
  #           In this case, the constant is obtained by extrapolating to the
  #           left from the `beta` percentile beyond the minimum value of `y`;
  #           the extrapolation is by a (small, positive) amount `gamma`.
  #
  stats <- quantile(y, c(alpha, 1/2, 1-alpha))
  code <- 0
  if (diff(diff(stats)==0)) {
    const <- Inf
    code <- 1
  }
  else {
    const <- (stats[2]^2 - stats[1]*stats[3]) / (stats[1] + stats[3] - 2*stats[2])
    y.min <- min(y)
    if (const < -y.min) {
      const <- gamma * quantile(y, beta) - (1+gamma) * y.min
      code <- 2
    }
  }
  list(offset=as.numeric(const), code=code)
}

#former f4
function_flexible <- function(par, focal_data, log_offset=0) {
	varepsilon_0 <- par[1]
	k <- par[2]
	a <- par[3]
	return(log(log_offset + varepsilon_0 + k*(focal_data$time)^a) - log(log_offset + focal_data$time))
}

generate_all_models <- function() {
	param_possibilities <- expand.grid(e=c(0, "e"), k=c(0, "k"), a=c(1, -1, "a")) #NA means it's free to vary
	param_possibilities$description <- paste0(
		ifelse(param_possibilities$a=="a", "free", ifelse(param_possibilities$a==-1, "hyperbola", "linear")),
		"_",
		param_possibilities$e,
		param_possibilities$k,
		ifelse(param_possibilities$a=="a", "a", "")
		)
	param_possibilities <- subset(param_possibilities, !(param_possibilities$description %in% c("linear_e0", "linear_00", "hyperbola_00", "free_00a", "free_e0a")))
	param_possibilities$nfreeparam <- apply(param_possibilities[,1:3], 1, function(x) {sum(is.na(suppressWarnings(as.numeric(x))))})
	return(param_possibilities)
}

optimize_rate_model<- function(focal_data, function_name, nparams, lb=-Inf, ub=Inf, nstarts_extra=10, all_algorithms=c("NLOPT_LN_BOBYQA", "NLOPT_LN_SBPLX", "NLOPT_LN_NEWUOA_BOUND"), log_offset=0) {
	model_distance_calls <<- 0
	model_distance_not_finite <<- 0
	par=rep(1, nparams)
	if(any(is.finite(ub))) {
		par[is.finite(ub)] <- ub[is.finite(ub)]
	}
	nfreeparams <- sum(!is.finite(ub))
	model_distance <- function(par, focal_data, log_offset) {
		model_distance_calls <<- model_distance_calls + 1
		predictions <- function_name(par, focal_data, log_offset)
		difference <- Inf
		difference <- sum((focal_data$log_rate_with_offset - predictions)^2)
		if(!is.finite(difference)) {
			model_distance_not_finite <<- model_distance_not_finite + 1
			#print(difference)
			difference <- 1e10
		}
		neglnL <- 0.5*nrow(focal_data)*log(difference/nrow(focal_data)) #yes, see lnL at https://en.wikipedia.org/wiki/Akaike_information_criterion#Comparison_with_least_squares, which is -0.5*n*log(RSS/n), so we get rid of the negative sign
		return(neglnL)
	}
	#return(optim(par=par, fn=model_distance, df=df, lower=lb, upper=ub, method="L-BFGS-B"))
	result <- nloptr::nloptr(x0=par, eval_f=model_distance,  lb=lb, ub=ub, opts=list(algorithm="NLOPT_LN_SBPLX", xtol_rel = 1e-4), focal_data=focal_data, log_offset=log_offset)
	
	# starting with lower param values, since they're often small
	par2 <- c(0.1, 0.0001, 1)[1:nparams]
	if(any(is.finite(ub))) {
		par2[is.finite(ub)] <- ub[is.finite(ub)]
	}
	
	result2 <- nloptr::nloptr(x0=par2, eval_f=model_distance,  lb=lb, ub=ub, opts=list(algorithm="NLOPT_LN_SBPLX", xtol_rel = 1e-4), focal_data=focal_data, log_offset=log_offset)
	
	if(result2$objective < result$objective) {
		result <- result2
	}
	
	for(start_index in sequence(nstarts_extra)) {
		par3 <- result$solution
		widths <- ub-lb
		sd_vector <- apply(rbind(widths, rep(0.1, length(widths))), 2, min)
		par3 <- rnorm(length(result$solution), mean=result$solution, sd=sd_vector)
		par3[par3<lb] <- lb[par3<lb]
		par3[par3>ub] <- ub[par3>ub]
		result3 <- nloptr::nloptr(x0=par3, eval_f=model_distance,  lb=lb, ub=ub, opts=list(algorithm=all_algorithms[1 + start_index%%length(all_algorithms)], xtol_rel = 1e-4), focal_data=focal_data, log_offset=log_offset)
		if(result3$objective < result$objective) {
			result <- result3
		}	
	}
	
	names(result$solution) <- c("e", "k", "a")[1:length(result$solution)]
	dentist_result <- dentist::dent_walk(par=result$solution, fn=model_distance, best_neglnL=result$objective, lower_bound=lb, upper_bound=ub, print_freq=1e6, focal_data=focal_data, log_offset=log_offset)
	result$dentist_result <- dentist_result
	result$nfreeparams <- nfreeparams
	print(paste0("model_distance_calls: ", model_distance_calls))
	print(paste0("model_distance_not_finite: ", model_distance_not_finite))
	return(result)
}

#' Check and clean input data
#' 
#' This will check and clean the input data. It expects column names of time, rate, and optionally citation; it can also include a numerator, denominator, and/or total_time columns. It will add missing columns.
#' @param all_data A data frame with columns of time, rate, and perhaps more.
#' @return A data frame with columns of time, rate, citation, numerator, denominator, and other columns in the input.
clean_input_data <- function(all_data) {
	if(!"citation" %in% colnames(all_data)) {
		all_data$citation <- "uncited"
	}
	if(!"log_offset" %in% colnames(all_data)) {
		all_data$log_offset <- 0
		for(dataset_name in unique(all_data$citation)) {
			indices <- which(all_data$citation==dataset_name)
			log_offset_result <- log.start(all_data$rate[indices])
			all_data$log_offset[indices] <- log_offset_result$offset	
		}
	}
	if(!"rate" %in% colnames(all_data)) {
		stop("Rate column not found in input data")
	}
	if(!"time" %in% colnames(all_data)) {
		stop("Time column not found in input data")
	}
	if(!"numerator" %in% colnames(all_data)) {
		all_data$numerator <- all_data$rate*all_data$time
	}
	if(!"denominator" %in% colnames(all_data)) {
		all_data$denominator <- all_data$time
	}
	if(!"total_time" %in% colnames(all_data)) {
		all_data$total_time <- NA
	}
	return(all_data)
}

#' Optimize rate model
#' 
#' This will dredge all the rate models for a given dataset. It expects column names of time, rate, and citation; optionally, it can also include a numerator, denominator, and/or total_time columns. It will return a list of results, one for each model and dataset.
#' 
#' By default it will use all models, but you can specify a subset of models to use.
#' @param all_data A data frame with columns of time, rate, and citation.
#' @param models A data frame with columns of e, k, a, and description
#' @param epsilon_lower The lower bound for the epsilon parameter.
#' @return A list of results, one for each model and dataset.
#' @export
optimization_over_all_data <- function(all_data, models=generate_all_models(), epsilon_lower=0) {
	all_data <- clean_input_data(all_data)
	datasets <- unique(all_data$citation)
	results <- list()
	for(dataset in datasets) {
		focal_data <- subset(all_data, citation==dataset)
		log_offset <- focal_data$log_offset[1]
		focal_data$log_rate_with_offset <- log(focal_data$rate + log_offset)
		
		for(model_index in sequence(nrow(models))) {
			lb <- c(epsilon_lower, rep(-Inf, 2))
			ub <- rep(Inf, 3)
			if(models$e[model_index]!="e") {
				lb[1] <- as.numeric(as.character(models$e[model_index]))
				ub[1] <- as.numeric(as.character(models$e[model_index]))
			}
			if(models$k[model_index]!="k") {
				lb[2] <- as.numeric(as.character(models$k[model_index]))
				ub[2] <- as.numeric(as.character(models$k[model_index]))
			}
			if(models$a[model_index]!="a") {
				lb[3] <- as.numeric(as.character(models$a[model_index]))
				ub[3] <- as.numeric(as.character(models$a[model_index]))
			}
			print(paste0("Optimizing model ", model_index, " of ", nrow(models), " ", models$description[model_index],  " for dataset ", dataset))
			print(paste0("lb: ", paste0(lb, collapse=", ")))
			print(paste0("ub: ", paste0(ub, collapse=", ")))
			
			local_result <- optimize_rate_model(focal_data, function_flexible, nparams=3, lb=lb, ub=ub, log_offset=log_offset)
			local_result$model <- models$description[model_index]
			local_result <- summarize_model(local_result, focal_data, function_flexible, log_offset=log_offset)
			results[[length(results)+1]] <- local_result
			names(results)[length(results)] <- paste0(dataset, "__", local_result$model)	

		}
	}	
	return(results)
}



summarize_model <- function(local_result, focal_data, function_name, log_offset) {
	local_result$n <- nrow(focal_data)
	local_result$AIC <- nrow(focal_data)*local_result$objective + 2*local_result$nfreeparams
	local_result$numerator <- focal_data$numerator
	local_result$denominator <- focal_data$denominator
	local_result$total_time <- focal_data$total_time
	local_result$time <- focal_data$time
	solution <- local_result$solution
	solution_nomserr <- solution
	solution_nomserr[1] <- 0
	local_result$predicted_log_rate_with_offset <- function_name(local_result$solution, focal_data, log_offset=log_offset)
	local_result$predicted_rate <- exp(local_result$predicted_log_rate_with_offset) - log_offset
	local_result$empirical_log_rate <- focal_data$log_rate
	local_result$empirical_log_rate_with_offset <- log(focal_data$rate + log_offset)
	local_result$empirical_rate <- focal_data$rate
	local_result$predicted_log_rate_with_offset_no_mserr <- rep(NA, length(local_result$empirical_log_rate))
	try({ local_result$predicted_log_rate_with_offset_no_mserr <- function_name(solution_nomserr, focal_data, log_offset=log_offset)})
	local_result$predicted_rate_no_mserr <- NA
	try({ local_result$predicted_rate_no_mserr <- exp(local_result$predicted_log_rate_with_offset_no_mserr) - log_offset})
	local_result$error_only_log_rate_with_offset <- local_result$predicted_log_rate_with_offset - local_result$predicted_log_rate_with_offset_no_mserr
	parameters_no_epsilon <- local_result$par
	#parameters_no_epsilon[1] <- 0
	#local_result$predicted_log_rate_no_mserr <- function_name(parameters_no_epsilon, df)
	return(local_result)
}

#' Generate summary information for all fitted models
#' 
#' This will generate a data frame with summary information for all fitted models.
#' @param minimization_approach_result A list of results from optimization_over_all_data
#' @return A data frame with summary information for all fitted models.
#' @export
summarize_all_fitted_models <- function(minimization_approach_result) {
	results <- data.frame()
	for (focal_model in names(minimization_approach_result)) {
		data_name <- strsplit(focal_model, "__")[[1]][1]
		focal_result <- minimization_approach_result[[focal_model]]
		params <- rep(NA,3)
		#solution <- focal_result$solution
		solution <- getElement(focal_result, "solution")
		params[1:length(solution)] <- solution
		names(params) <- c("e", "k", "a")
		if(ncol(focal_result$dentist_result$all_ranges)<3) { # pad to handle only getting 1 or 2 params
			focal_result$dentist_result$all_ranges <- cbind(focal_result$dentist_result$all_ranges, rep(NA, nrow(focal_result$dentist_result$all_ranges)))
			focal_result$dentist_result$all_ranges <- cbind(focal_result$dentist_result$all_ranges, rep(NA, nrow(focal_result$dentist_result$all_ranges)))
		}
		focal_df <- suppressWarnings(data.frame(
			dataset=data_name, 
			model=focal_result$model, 
			n=focal_result$n, 
			AIC=focal_result$AIC, 
			objective=focal_result$objective, 
			nfreeparams=focal_result$nfreeparams, 
			param_e=params['e'], 
			param_k=params['k'], 
			param_a=params['a'], 
			param_e_lower = focal_result$dentist_result$all_ranges['lower.CI', 1], 
			param_e_upper =  focal_result$dentist_result$all_ranges['upper.CI', 1], 
			param_k_lower = focal_result$dentist_result$all_ranges['lower.CI', 2], 
			param_k_upper =  focal_result$dentist_result$all_ranges['upper.CI', 2], 
			param_a_lower = focal_result$dentist_result$all_ranges['lower.CI', 3], 
			param_a_upper =  focal_result$dentist_result$all_ranges['upper.CI', 3], predicted_log_rate_with_offset=focal_result$predicted_log_rate_with_offset, 
			empirical_log_rate=focal_result$empirical_log_rate, 
			empirical_log_rate_with_offset=focal_result$empirical_log_rate_with_offset,
			predicted_log_rate_with_offset_no_mserr=focal_result$predicted_log_rate_with_offset_no_mserr, 
			error_only_log_rate_with_offset=focal_result$error_only_log_rate_with_offset, 
			predicted_rate=focal_result$predicted_rate,
			empirical_rate=focal_result$empirical_rate,
			predicted_rate_no_mserr=focal_result$predicted_rate_no_mserr,
			time=focal_result$time, 
			numerator=focal_result$numerator, 
			total_time=focal_result$total_time, 
			denominator=focal_result$denominator
		))
		focal_df_tall <- focal_df |> tidyr::pivot_longer(cols=c("predicted_log_rate_with_offset", "empirical_log_rate_with_offset", "predicted_log_rate_with_offset_no_mserr"), names_to="log_rate_type", values_to="log_rate") |> tidyr::pivot_longer(cols=c("predicted_rate", "empirical_rate", "predicted_rate_no_mserr"), names_to="rate_type", values_to="rate")

		#focal_df_tall <- focal_df_tall |> tidyr::pivot_longer(cols=c("predicted_nonlog_rate", "empirical_nonlog_rate", "predicted_nonlog_rate_no_mserr", "error_only_nonlog_rate"), names_to="rate_type", values_to="nonlog_rate")
		results <- rbind(results, focal_df_tall)
	}
	results$deltaAIC <- NA
	for (focal_dataset in unique(results$dataset)) {
		focal_rows <- which(results$dataset==focal_dataset)
		results$deltaAIC[focal_rows] <- results$AIC[focal_rows] - min(results$AIC[focal_rows])
	}
	return(results)	
}

summarize_all_models <- function(minimization_approach_result_summarized) {
	models <- minimization_approach_result_summarized |> dplyr::select(-rate_type) |> dplyr::select(-log_rate) |> dplyr::select(-log_time) |> dplyr::select(-nparams) |> dplyr::distinct()
	models <- models[order(models$dataset, models$deltaAIC),]
	return(models)
}


optimization_and_summarization_over_randomized_data <- function(all_data, nreps=5, epsilon_lower=0) {
	final_result <- data.frame()
	for(rep_index in sequence(nreps)) {
		local_result <- summarize_all_fitted_models(optimization_over_all_data(randomize_within_dataset(all_data), epsilon_lower=epsilon_lower))
		local_result$rep <- rep_index
		final_result <- rbind(final_result, local_result)	
	}
	return(final_result)
}

merge_original_and_random <- function(minimization_models_summarized, randomized_data_models_summarized) {
	minimization_models_summarized$rep <- "Original"
	randomized_data_models_summarized$rep <- paste0("Rep ", randomized_data_models_summarized$rep)
	merged <- dplyr::bind_rows(minimization_models_summarized, randomized_data_models_summarized)
	return(merged)
}

merge_best_original_and_random <- function(minimization_models_best, randomized_data_models_best) {
	minimization_models_best$rep <- "Original"
	randomized_data_models_best$rep <- paste0("Rep ", randomized_data_models_best$rep)
	merged <- dplyr::bind_rows(minimization_models_best, randomized_data_models_best)
	return(merged)
}

randomize_within_dataset <- function(all_data) {
	all_data <- clean_input_data(all_data)
	result <- data.frame()
	datasets <- unique(all_data$citation)
	for (dataset in datasets) {
		print(dataset)
		focal_data <- subset(all_data, citation==dataset)
		focal_data$rate <- sample(focal_data$numerator)/focal_data$denominator
		result <- rbind(result, focal_data)
	}
	return(result)
}
