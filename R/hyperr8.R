#' Run all analyses
#' 
#' By default will run analyses on all datasets. It expects an input data frame with columns of time and rate, each row representing a single observation. To handle multiple datasets at once, have a citation column indicating the name for the dataset (which could be a full citation or as simple as "dataset1", "dataset2", etc.). You can also include a column for numerator (the number of events) and denominator (the number of opportunities for events). If you have a column for total_time, it will use that for randomization rather than time (this is important for things like phylogenetic trees, where the time over which events happen is not just the overall age of the tree but the sum of its branch lengths). 
#' 
#' Data can be randomized to help suggest whether an observed pattern is spurious. 
#' @param all_data A data frame with columns of time, rate, and perhaps citation, numerator, denominator, and/or total_time (other columns will be ignored).
#' @param nreps The number of times to randomize the data within each dataset.
#' @param epsilon_lower The lower bound for the epsilon parameter.
#' @param nstep_dentist The number of steps for the dentist algorithm.
#' @return A data.frame of results with class hyperr8.
#' @export
#' @examples
#' library(hyperr8)
#' car_data <- generate_car_simulation()
#' all_run <- hyperr8_run(car_data)
#' library(ggplot2)
#' plot(all_run)
hyperr8_run <- function(all_data, nreps=5, epsilon_lower=-Inf, nstep_dentist=1000) {
	all_data <- clean_input_data(all_data)
	original <- summarize_all_fitted_models(optimization_over_all_data(all_data, epsilon_lower=epsilon_lower, nstep_dentist=nstep_dentist))
	original$rep <- "Original"
	if(nreps>0) {
		randomized <- optimization_and_summarization_over_randomized_data(all_data, nreps=nreps, epsilon_lower=epsilon_lower, nstep_dentist=nstep_dentist)
		randomized$rep <- paste0("Rep ", randomized$rep)
		merged <- dplyr::bind_rows(original, randomized)
	} else {
		merged <- original
	}
	merged <- as.data.frame(merged) # I am a Klingon with bad spelling. Death to tibbles. They have no honor.
	class(merged) <- c("hyperr8", class(merged))
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
plot.hyperr8 <- function(x, loglog=TRUE, ...) {
	x <- x |> tidyr::pivot_longer(cols=c("predicted_log_rate_with_offset", "empirical_log_rate_with_offset", "predicted_log_rate_with_offset_no_mserr"), names_to="log_rate_type", values_to="log_rate") |> tidyr::pivot_longer(cols=c("predicted_rate", "empirical_rate", "predicted_rate_no_mserr"), names_to="rate_type", values_to="rate")
	npoints <- nrow(subset(x, rate_type=='empirical_rate' & deltaAIC==0))
	#alpha <- max(0.005, min(0.7, 10/sqrt(npoints)))
	alpha <- 0.1
	x$data_rep <- paste0(x$dataset, "\n", x$rep)
	gcool <- ggplot(subset(x, rate_type=='empirical_rate' & deltaAIC==0), aes(x=time, y=rate)) + geom_point(alpha=alpha, shape=20) + facet_wrap(~data_rep, scales="free") + theme_bw() + xlab("Time") + ylab("Rate")
	if(loglog) {
		gcool <- gcool + scale_x_continuous(trans = "log") + scale_y_continuous(trans = "log")
	}
	gcool <- gcool + geom_line(data=subset(x, rate_type=='predicted_rate' & deltaAIC==0), aes(x=time, y=rate, group=rep, colour=model))
	gcool
	return(gcool)
}

#' Summarize hyperr8 results
#' 
#' This will summarize the results from hyperr8_run.
#' @param object The output from hyperr8_run.
#' @param ... Other arguments to pass to the summary function.
#' @return A list with summary information.
#' @export
summary.hyperr8 <- function(object, ...) {
	class(object) <- 'data.frame'
	distinct_df <- dplyr::distinct(object, dataset, model, n, objective, nfreeparams, param_h, param_m, param_b, param_h_lower, param_h_upper, param_m_lower, param_m_upper, param_b_lower, param_b_upper, deltaAIC, rep)
	original <- subset(distinct_df, rep=="Original")
	original_best <- subset(original, deltaAIC==0)
	return_list <- list(original=original, original_best=original_best)
	randomized <- subset(distinct_df, rep!="Original")
	if(nrow(randomized)>0) {
		randomized_best <- subset(randomized, deltaAIC==0)
		return_list$randomized <- randomized
		return_list$randomized_best <- randomized_best
	}
	return(return_list)
}

#' Print hyperr8 results
#' 	
#' This will summarize the results from hyperr8_run.
#' @param x The output from hyperr8_run.
#' @param ... Other arguments to pass to the printing function.
#' @return A data.frame with summary information.
#' @export
print.hyperr8 <- function(x,...) {
	distinct_df <- dplyr::distinct(x, dataset, model, n, objective, nfreeparams, param_h, param_m, param_b, param_h_lower, param_h_upper, param_m_lower, param_m_upper, param_b_lower, param_b_upper, deltaAIC, rep)
	original <- subset(distinct_df, rep=="Original")
	original <- original[order(original$deltaAIC),]
	return(as.data.frame(original))
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


# rate = e/t + m*t^a + b

# This returns the log rate.
function_flexible <- function(par, focal_data, log_offset=0, do_log=TRUE, do_dnorm=FALSE) {
	h <- par[1]
	m <- par[2]
	b <- par[3]
	time <- focal_data$time
	result <- log_offset +
		h/time +
		m*time +
		b
	if(do_dnorm) {
		result <- log_offset + 
			(h*sqrt(2/pi) + # folded normal
			m*time^2 +
			b*time) / time
			
	}
	if(do_log) {
		result <- log(result)
	}
	return(result)
}

# rate = e/t + m*t^a + b

generate_all_models <- function() {
	param_possibilities <- expand.grid(h=c(0, "h"), m=c(0, "m"), b=c(0, "b")) 
	param_possibilities$description <- paste0(
		param_possibilities$h,
		param_possibilities$m,
		param_possibilities$b
	)
	
	param_possibilities <- subset(param_possibilities, !(param_possibilities$description %in% c("000"))) # rates are 0

	
	

	rownames(param_possibilities) <- NULL
	
	
	
	param_possibilities$nfreeparam <- apply(param_possibilities[,1:3], 1, function(x) {sum(is.na(suppressWarnings(as.numeric(x))))})
	return(param_possibilities)
}

generate_all_models_OLD <- function() {
	param_possibilities <- expand.grid(e=c(0, "e"), k=c(0, "k"), a=c(1, -1)) #NA means it's free to vary
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

optimize_rate_model<- function(focal_data, function_name, lb=-Inf, ub=Inf, nstarts_extra=10, all_algorithms=c("NLOPT_LN_BOBYQA", "NLOPT_LN_SBPLX", "NLOPT_LN_NEWUOA_BOUND"), log_offset=0, nstep_dentist=1000) {
	model_distance_calls <<- 0
	model_distance_not_finite <<- 0
	par=c(rep(.1, 2), mean(focal_data$rate))
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
	par2 <- c(0.01, 0.0001, 0.1*mean(focal_data$rate))
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
		par3 <- stats::rnorm(length(result$solution), mean=result$solution, sd=sd_vector)
		par3[par3<lb] <- lb[par3<lb]
		par3[par3>ub] <- ub[par3>ub]
		result3 <- nloptr::nloptr(x0=par3, eval_f=model_distance,  lb=lb, ub=ub, opts=list(algorithm=all_algorithms[1 + start_index%%length(all_algorithms)], xtol_rel = 1e-4), focal_data=focal_data, log_offset=log_offset)
		if(result3$objective < result$objective) {
			result <- result3
		}	
	}
	
	names(result$solution) <- c("h", "m", "b")[1:length(result$solution)]
	dentist_result <- suppressWarnings({dentist::dent_walk(par=result$solution, fn=model_distance, best_neglnL=result$objective, lower_bound=lb, upper_bound=ub, print_freq=1e6, focal_data=focal_data, log_offset=log_offset, nsteps=nstep_dentist)}) # if the likelihood improves during dent_walk, it will return the improved solution
	if(min(dentist_result$results$neglnL) < result$objective) {
		best_row <- which.min(dentist_result$results$neglnL)
		result$solution <- unname(unlist(dentist_result$results[best_row,-1]))
		result$objective <-  unname(unlist(dentist_result$results[best_row,1]))
	}
	result$dentist_result <- dentist_result
	result$nfreeparams <- nfreeparams
	result$datum_id <- focal_data$datum_id
	#(paste0("model_distance_calls: ", model_distance_calls))
	#print(paste0("model_distance_not_finite: ", model_distance_not_finite))
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
	if(!"datum_id" %in% colnames(all_data)) {
		all_data$datum_id <- sequence(nrow(all_data))
	}
	if(!"log_offset" %in% colnames(all_data)) {
		all_data$log_offset <- 0
		for(dataset_name in unique(all_data$citation)) {
			indices <- which(all_data$citation==dataset_name)
			if(any(all_data$rate[indices]<0)) {
				stop("Negative rates found in input data")
			}
			if(any(all_data$rate[indices]==0)) {
				all_data$log_offset[indices] <- 1 # log1p	
			}
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
#' @param nstep_dentist The number of steps for the dentist algorithm.
#' @return A list of results, one for each model and dataset.
#' @export
optimization_over_all_data <- function(all_data, models=generate_all_models(), epsilon_lower=-Inf, nstep_dentist=1000) {
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
			if(models$h[model_index]!="h") {
				lb[1] <- as.numeric(as.character(models$h[model_index]))
				ub[1] <- as.numeric(as.character(models$h[model_index]))
			}
			if(models$m[model_index]!="m") {
				lb[2] <- as.numeric(as.character(models$m[model_index]))
				ub[2] <- as.numeric(as.character(models$m[model_index]))
			}
			if(models$b[model_index]!="b") {
				lb[3] <- as.numeric(as.character(models$b[model_index]))
				ub[3] <- as.numeric(as.character(models$b[model_index]))
			}
			print(paste0("Optimizing model ", model_index, " of ", nrow(models), " ", models$description[model_index],  " for dataset ", dataset))
		#	print(paste0("lb: ", paste0(lb, collapse=", ")))
		#	print(paste0("ub: ", paste0(ub, collapse=", ")))
			
			local_result <- optimize_rate_model(focal_data, function_flexible, lb=lb, ub=ub, log_offset=log_offset, nstep_dentist=nstep_dentist)
			local_result$model <-  models$description[model_index]
			local_result$description <- models$description[model_index]
			local_result <- summarize_model(local_result, focal_data, function_flexible, log_offset=log_offset)
			results[[length(results)+1]] <- local_result
			names(results)[length(results)] <- paste0(dataset, "__", local_result$model)	

		}
	}	
	return(results)
}



summarize_model <- function(local_result, focal_data, function_name, log_offset, do_dnorm=FALSE) {
	local_result$n <- nrow(focal_data)
	local_result$AIC <- 2*local_result$objective + 2*local_result$nfreeparams
	local_result$numerator <- focal_data$numerator
	local_result$denominator <- focal_data$denominator
	local_result$total_time <- focal_data$total_time
	local_result$time <- focal_data$time
	solution <- local_result$solution
	solution_nomserr <- solution
	solution_nomserr[1] <- 0
	local_result$predicted_log_rate_with_offset <- function_name(local_result$solution, focal_data, log_offset=log_offset, do_dnorm=do_dnorm)
	local_result$predicted_rate <- exp(local_result$predicted_log_rate_with_offset) - log_offset
	local_result$empirical_log_rate <- focal_data$log_rate
	local_result$empirical_log_rate_with_offset <- log(focal_data$rate + log_offset)
	local_result$empirical_rate <- focal_data$rate
	local_result$hyperbolic_component <- solution[1]/focal_data$time
	local_result$linear_component <- solution[2]*focal_data$time
	local_result$constant_component <- solution[3]+0*focal_data$time
	
	local_result$hyperbolic_component_proportion <- abs(local_result$hyperbolic_component)/(abs(local_result$hyperbolic_component) + abs(local_result$linear_component) + abs(local_result$constant_component))
	local_result$linear_component_proportion <- abs(local_result$linear_component)/(abs(local_result$hyperbolic_component) + abs(local_result$linear_component) + abs(local_result$constant_component))
	local_result$constant_component_proportion <- abs(local_result$constant_component)/(abs(local_result$hyperbolic_component) + abs(local_result$linear_component) + abs(local_result$constant_component))
	
	
	local_result$predicted_log_rate_with_offset_no_mserr <- NA
	local_result$predicted_rate_with_offset_no_mserr <- NA
	local_result$predicted_rate_no_mserr <- NA

	try({ local_result$predicted_rate_with_offset_no_mserr <- function_name(solution_nomserr, focal_data, log_offset=log_offset, do_log=FALSE)})
	try({ local_result$predicted_rate_no_mserr <- local_result$predicted_rate_with_offset_no_mserr - log_offset})
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
		names(params) <- c("h", "m", "b")
		while(ncol(focal_result$dentist_result$all_ranges)<3) { # pad to handle only getting 1:2 params
			focal_result$dentist_result$all_ranges <- cbind(focal_result$dentist_result$all_ranges, rep(NA, nrow(focal_result$dentist_result$all_ranges)))
		}
		focal_df <- suppressWarnings(data.frame(
			dataset=data_name, 
			model=focal_result$model, 
			n=focal_result$n, 
			AIC=focal_result$AIC, 
			objective=focal_result$objective, 
			nfreeparams=focal_result$nfreeparams, 
			param_h=params['h'], 
			param_m=params['m'], 
			param_b=params['b'],
			param_h_lower = focal_result$dentist_result$all_ranges['lower.CI', 1], 
			param_h_upper =  focal_result$dentist_result$all_ranges['upper.CI', 1], 
			param_m_lower = focal_result$dentist_result$all_ranges['lower.CI', 2], 
			param_m_upper =  focal_result$dentist_result$all_ranges['upper.CI', 2], 
			param_b_lower = focal_result$dentist_result$all_ranges['lower.CI', 3],
			param_b_upper =  focal_result$dentist_result$all_ranges['upper.CI', 3],
			predicted_log_rate_with_offset=focal_result$predicted_log_rate_with_offset, 
			empirical_log_rate=focal_result$empirical_log_rate, 
			empirical_log_rate_with_offset=focal_result$empirical_log_rate_with_offset,
			predicted_log_rate_with_offset_no_mserr=focal_result$predicted_log_rate_with_offset_no_mserr, 
			predicted_rate=focal_result$predicted_rate,
			empirical_rate=focal_result$empirical_rate,
			predicted_rate_no_mserr=focal_result$predicted_rate_no_mserr,
			time=focal_result$time, 
			hyperbolic_component=focal_result$hyperbolic_component,
			linear_component=focal_result$linear_component,
			constant_component=focal_result$constant_component,
			hyperbolic_component_proportion=focal_result$hyperbolic_component_proportion,
			linear_component_proportion=focal_result$linear_component_proportion,
			constant_component_proportion=focal_result$constant_component_proportion,
			numerator=focal_result$numerator, 
			total_time=focal_result$total_time, 
			denominator=focal_result$denominator,
			datum_id=focal_result$datum_id
		))
	#	focal_df_tall <- focal_df |> tidyr::pivot_longer(cols=c("predicted_log_rate_with_offset", "empirical_log_rate_with_offset", "predicted_log_rate_with_offset_no_mserr"), names_to="log_rate_type", values_to="log_rate") |> tidyr::pivot_longer(cols=c("predicted_rate", "empirical_rate", "predicted_rate_no_mserr"), names_to="rate_type", values_to="rate")
		
	#	focal_df_tall <- focal_df_tall |> tidyr::pivot_longer(cols=c("hyperbolic_component", "linear_component", "constant_component"), names_to="component_type", values_to="component")
	#	focal_df_tall$component_type <- gsub("_component", "", focal_df_tall$component_type)
		
	#	focal_df_tall <- focal_df_tall |> tidyr::pivot_longer(cols=c("hyperbolic_component_proportion", "linear_component_proportion", "constant_component_proportion"), names_to="component_proportion_type", values_to="component_proportion")
	#	focal_df_tall$component_proportion_type <- gsub("_component_proportion", "", focal_df_tall$component_proportion_type)

	#	results <- rbind(results, focal_df_tall)
		results <- rbind(results, focal_df)
	}
	results$deltaAIC <- NA
	for (focal_dataset in unique(results$dataset)) {
		focal_rows <- which(results$dataset==focal_dataset)
		results$deltaAIC[focal_rows] <- results$AIC[focal_rows] - min(results$AIC[focal_rows])
	}
	return(results)	
}

summarize_all_models <- function(minimization_approach_result_summarized) {
	models <- minimization_approach_result_summarized |> dplyr::select(-rate_type) |> dplyr::select(-log_rate) |> dplyr::select(-log_time) |> dplyr::select(-3) |> dplyr::distinct()
	models <- models[order(models$dataset, models$deltaAIC),]
	return(models)
}

#' Plot contributions of each component
#' 
#' This will plot the contributions of each component.
#' @param x The output from summarize_all_fitted_models, filtered to a single dataset and model
#' @param scaling_factor The scaling factor for the plot, which will be used to determine the height of the ribbons.
#' @return A ggplot2 object.
#' @export
plot_proportion_with_offset <- function(x, scaling_factor=0.3) {
	x <- x |> tidyr::pivot_wider(names_from="component_proportion_type", values_from="component_proportion")
	x <- subset(x, rate_type=="predicted_rate")
	starting_values <- log(x$rate)-0.5*scaling_factor
	x$constant_component_proportion_ribbon_lower <- starting_values
	x$constant_component_proportion_ribbon_upper <- x$constant_component_proportion_ribbon_lower+scaling_factor*(x$constant)
	x$linear_component_proportion_ribbon_lower <- x$constant_component_proportion_ribbon_upper
	x$linear_component_proportion_ribbon_upper <- x$linear_component_proportion_ribbon_lower + scaling_factor*(x$linear)
	x$hyperbolic_component_proportion_ribbon_lower <- x$linear_component_proportion_ribbon_upper 
	x$hyperbolic_component_proportion_ribbon_upper <- x$hyperbolic_component_proportion_ribbon_lower + scaling_factor*(x$hyperbolic)
	
	x <- dplyr::distinct(x, dataset, time, rate, rate_type, constant_component_proportion_ribbon_lower, constant_component_proportion_ribbon_upper, linear_component_proportion_ribbon_lower, linear_component_proportion_ribbon_upper, hyperbolic_component_proportion_ribbon_lower, hyperbolic_component_proportion_ribbon_upper)
	
	# now do geom_ribbon
	gcool <- ggplot2::ggplot(x, aes(x=time, y=rate)) + 
		geom_ribbon(aes(ymin=constant_component_proportion_ribbon_lower, ymax=constant_component_proportion_ribbon_upper), fill="red") + 
		geom_ribbon(aes(ymin=linear_component_proportion_ribbon_lower, ymax=linear_component_proportion_ribbon_upper), fill="blue") + 
		geom_ribbon(aes(ymin=hyperbolic_component_proportion_ribbon_lower, ymax=hyperbolic_component_proportion_ribbon_upper), fill="green") 
		#gcool <- gcool + geom_line(aes(y=log(rate)),colour="black") 
		gcool <- gcool + scale_x_continuous(trans="log")
		gcool <- gcool + scale_y_continuous(labels=function(y) exp(y))
	
	return(gcool)
}

optimization_and_summarization_over_randomized_data <- function(all_data, nreps=5, epsilon_lower=-Inf, nstep_dentist=1000) {
	final_result <- data.frame()
	for(rep_index in sequence(nreps)) {
		local_result <- summarize_all_fitted_models(optimization_over_all_data(randomize_within_dataset(all_data), epsilon_lower=epsilon_lower, nstep_dentist=nstep_dentist))
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

#' Run all analyses using norm appoach
#' 
#' By default will run analyses on all datasets. It expects an input data frame with columns of time and rate, each row representing a single observation. To handle multiple datasets at once, have a citation column indicating the name for the dataset (which could be a full citation or as simple as "dataset1", "dataset2", etc.). You can also include a column for numerator (the number of events) and denominator (the number of opportunities for events). If you have a column for total_time, it will use that for randomization rather than time (this is important for things like phylogenetic trees, where the time over which events happen is not just the overall age of the tree but the sum of its branch lengths). 
#' 
#' Data can be randomized to help suggest whether an observed pattern is spurious. 
#' @param all_data A data frame with columns of time, rate, and perhaps citation, numerator, denominator, and/or total_time (other columns will be ignored).
#' @param nreps The number of times to randomize the data within each dataset.
#' @param epsilon_lower The lower bound for the epsilon parameter.
#' @param nstep_dentist The number of steps for the dentist algorithm.
#' @return A data.frame of results with class hyperr8.
#' @export
#' @examples
#' library(hyperr8); car_data <- generate_car_simulation(); all_run <- hyperr8_norm_run(car_data)
#' subset(as.data.frame(all_run), datum_id==1)
#' library(ggplot2)
#' plot(all_run)
hyperr8_norm_run <- function(all_data, nreps=5, epsilon_lower=-Inf, nstep_dentist=1000) {
	all_data <- clean_input_data(all_data)
	original <- summarize_all_fitted_models_norm_approach(optimization_norm_over_all_data(all_data, epsilon_lower=epsilon_lower, nstep_dentist=nstep_dentist))
	original$rep <- "Original"
	# if(nreps>0) {
	# 	randomized <- optimization_and_summarization_over_randomized_data(all_data, nreps=nreps, epsilon_lower=epsilon_lower, nstep_dentist=nstep_dentist)
	# 	randomized$rep <- paste0("Rep ", randomized$rep)
	# 	merged <- dplyr::bind_rows(original, randomized)
	# } else {
		merged <- original
	# }
	merged <- as.data.frame(merged) # I am a Klingon with bad spelling. Death to tibbles. They have no honor.
	class(merged) <- c("hyperr8", class(merged))
	return(merged)
}

#' Optimize rate model using noise model
#' 
#' This will dredge all the rate models for a given dataset. It expects column names of time, rate, and citation; optionally, it can also include a numerator, denominator, and/or total_time columns. It will return a list of results, one for each model and dataset.
#' 
#' By default it will use all models, but you can specify a subset of models to use.
#' @param all_data A data frame with columns of time, rate, and citation.
#' @param models A data frame with columns of e, k, a, and description
#' @param epsilon_lower The lower bound for the epsilon parameter.
#' @param nstep_dentist The number of steps for the dentist algorithm.
#' @return A list of results, one for each model and dataset.
#' @export
optimization_norm_over_all_data <- function(all_data, models=generate_all_models()[c(3,5,7),], epsilon_lower=-Inf, nstep_dentist=1000) {
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
			if(models$h[model_index]!="h") {
				lb[1] <- as.numeric(as.character(models$h[model_index]))
				ub[1] <- as.numeric(as.character(models$h[model_index]))
			}
			if(models$m[model_index]!="m") {
				lb[2] <- as.numeric(as.character(models$m[model_index]))
				ub[2] <- as.numeric(as.character(models$m[model_index]))
			}
			if(models$b[model_index]!="b") {
				lb[3] <- as.numeric(as.character(models$b[model_index]))
				ub[3] <- as.numeric(as.character(models$b[model_index]))
			}
			print(paste0("Optimizing model ", model_index, " of ", nrow(models), " ", models$description[model_index],  " for dataset ", dataset))
		#	print(paste0("lb: ", paste0(lb, collapse=", ")))
		#	print(paste0("ub: ", paste0(ub, collapse=", ")))
			
			local_result <- optimize_norm_rate_model(focal_data, function_flexible, lb=lb, ub=ub, log_offset=log_offset, nstep_dentist=nstep_dentist)
			local_result$model <-  models$description[model_index]
			local_result$description <- models$description[model_index]
			local_result <- summarize_model(local_result, focal_data, function_flexible, log_offset=log_offset, do_dnorm=TRUE)
			results[[length(results)+1]] <- local_result
			names(results)[length(results)] <- paste0(dataset, "__", local_result$model)	
			print(local_result)

		}
	}	
	return(results)
}



optimize_norm_rate_model<- function(focal_data, function_name, lb=-Inf, ub=Inf, nstarts_extra=5, all_algorithms=c("NLOPT_LN_BOBYQA", "NLOPT_LN_SBPLX", "NLOPT_LN_NEWUOA_BOUND"), log_offset=0, nstep_dentist=1000) {
	model_distance_calls <<- 0
	model_distance_not_finite <<- 0
	par=c(100, 0.2, max(0.001,min(focal_data$rate)))
	if(any(is.finite(ub))) {
		par[is.finite(ub)] <- ub[is.finite(ub)]
	}
	nfreeparams <- sum(!is.finite(ub))
	

	likelihood_function <- function(par, focal_data) {
		par_no_h <- par
		par_no_h[1] <- 0
		means <- function_name(par_no_h, focal_data, 0, do_log=FALSE)
		negLnL <- -1 * sum(stats::dnorm(focal_data$rate*focal_data$time, mean=means*focal_data$time, sd=par[1], log=TRUE))
		return(negLnL)
	}
	
	#return(optim(par=par, fn=model_distance, df=df, lower=lb, upper=ub, method="L-BFGS-B"))
	result <- nloptr::nloptr(x0=par, eval_f=likelihood_function,  lb=lb, ub=ub, opts=list(algorithm="NLOPT_LN_SBPLX", xtol_rel = 1e-4), focal_data=focal_data)
	
	
	if(1==1) {
	# starting with lower param values, since they're often small
		# par2 <- c(result$solution*runif(length(result$solution), 0.9, 2))
		# if(any(is.finite(ub))) {
		# 	par2[is.finite(ub)] <- ub[is.finite(ub)]
		# }
		
		# result2 <- nloptr::nloptr(x0=par2, eval_f=likelihood_function,  lb=lb, ub=ub, opts=list(algorithm="NLOPT_LN_SBPLX", xtol_rel = 1e-4), focal_data=focal_data)
		
		# if(result2$objective < result$objective) {
		# 	result <- result2
		# }
		
		for(start_index in sequence(nstarts_extra)) {
			par3 <- c(result$solution*runif(length(result$solution), 0.9, 2))
			widths <- ub-lb
			print(paste(start_index, " ", all_algorithms[1 + start_index%%length(all_algorithms)]))

			#sd_vector <- apply(rbind(widths, rep(0.1, length(widths))), 2, min)
			#par3 <- stats::rnorm(length(result$solution), mean=result$solution, sd=sd_vector)
			par3[par3<lb] <- lb[par3<lb]
			par3[par3>ub] <- ub[par3>ub]
			result3 <- nloptr::nloptr(x0=par3, eval_f=likelihood_function,  lb=lb, ub=ub, opts=list(algorithm=all_algorithms[1 + start_index%%length(all_algorithms)], xtol_rel = 1e-4, maxtime=120, maxeval=100), focal_data=focal_data)
			if(result3$objective < result$objective) {
				result <- result3
			}	
			
			print(result3$objective)

		}
	}
	
	names(result$solution) <- c("h", "m", "b")[1:length(result$solution)]
	#dentist_result <- suppressWarnings({dentist::dent_walk(par=result$solution, fn=likelihood_function, best_neglnL=result$objective, lower_bound=lb, upper_bound=ub, print_freq=1e6, focal_data=focal_data, nsteps=nstep_dentist)}) # if the likelihood improves during dent_walk, it will return the improved solution
	# if(min(dentist_result$results$neglnL) < result$objective) {
	# 	best_row <- which.min(dentist_result$results$neglnL)
	# 	result$solution <- unname(unlist(dentist_result$results[best_row,-1]))
	# 	result$objective <-  unname(unlist(dentist_result$results[best_row,1]))
	# }
	# result$dentist_result <- dentist_result
	result$nfreeparams <- nfreeparams
	result$datum_id <- focal_data$datum_id
	#(paste0("model_distance_calls: ", model_distance_calls))
	#print(paste0("model_distance_not_finite: ", model_distance_not_finite))
	return(result)
}


#' Generate summary information for all fitted models
#' 
#' This will generate a data frame with summary information for all fitted models.
#' @param minimization_approach_result A list of results from optimization_over_all_data
#' @return A data frame with summary information for all fitted models.
#' @export
summarize_all_fitted_models_norm_approach <- function(minimization_approach_result) {
	results <- data.frame()
	for (focal_model in names(minimization_approach_result)) {
		data_name <- strsplit(focal_model, "__")[[1]][1]
		focal_result <- minimization_approach_result[[focal_model]]
		params <- rep(NA,3)
		#solution <- focal_result$solution
		solution <- getElement(focal_result, "solution")
		params[1:length(solution)] <- solution
		names(params) <- c("h", "m", "b")

		focal_df <- suppressWarnings(data.frame(
			dataset=data_name, 
			model=focal_result$model, 
			n=focal_result$n, 
			AIC=focal_result$AIC, 
			objective=focal_result$objective, 
			nfreeparams=focal_result$nfreeparams, 
			param_h=params['h'], 
			param_m=params['m'], 
			param_b=params['b'],
			predicted_log_rate_with_offset=focal_result$predicted_log_rate_with_offset, 
			empirical_log_rate=focal_result$empirical_log_rate, 
			empirical_log_rate_with_offset=focal_result$empirical_log_rate_with_offset,
			predicted_log_rate_with_offset_no_mserr=focal_result$predicted_log_rate_with_offset_no_mserr, 
			predicted_rate=focal_result$predicted_rate,
			empirical_rate=focal_result$empirical_rate,
			predicted_rate_no_mserr=focal_result$predicted_rate_no_mserr,
			time=focal_result$time, 
			hyperbolic_component=focal_result$hyperbolic_component,
			linear_component=focal_result$linear_component,
			constant_component=focal_result$constant_component,
			hyperbolic_component_proportion=focal_result$hyperbolic_component_proportion,
			linear_component_proportion=focal_result$linear_component_proportion,
			constant_component_proportion=focal_result$constant_component_proportion,
			numerator=focal_result$numerator, 
			total_time=focal_result$total_time, 
			denominator=focal_result$denominator,
			datum_id=focal_result$datum_id
		))
	#	focal_df_tall <- focal_df |> tidyr::pivot_longer(cols=c("predicted_log_rate_with_offset", "empirical_log_rate_with_offset", "predicted_log_rate_with_offset_no_mserr"), names_to="log_rate_type", values_to="log_rate") |> tidyr::pivot_longer(cols=c("predicted_rate", "empirical_rate", "predicted_rate_no_mserr"), names_to="rate_type", values_to="rate")
		
	#	focal_df_tall <- focal_df_tall |> tidyr::pivot_longer(cols=c("hyperbolic_component", "linear_component", "constant_component"), names_to="component_type", values_to="component")
	#	focal_df_tall$component_type <- gsub("_component", "", focal_df_tall$component_type)
		
	#	focal_df_tall <- focal_df_tall |> tidyr::pivot_longer(cols=c("hyperbolic_component_proportion", "linear_component_proportion", "constant_component_proportion"), names_to="component_proportion_type", values_to="component_proportion")
	#	focal_df_tall$component_proportion_type <- gsub("_component_proportion", "", focal_df_tall$component_proportion_type)

	#	results <- rbind(results, focal_df_tall)
		results <- rbind(results, focal_df)
	}
	results$deltaAIC <- NA
	for (focal_dataset in unique(results$dataset)) {
		focal_rows <- which(results$dataset==focal_dataset)
		results$deltaAIC[focal_rows] <- results$AIC[focal_rows] - min(results$AIC[focal_rows])
	}
	return(results)	
}


#' This is data to be included in the package
#'
#' @name yule_sim
#' @docType data
#' @keywords data
#' @usage data(yule_sim)
#' @format a data.frame with 25,000 rows and 8 columns
#' 
#' This has a dataset of estimated speciation rate for every tree simulated from a pure birth simulation with a speciation rate 0.1 and total time pulled from a distribution uniform on a log space ranging from 1 to 50 time units. Trees started with two taxa; trees with no speciation events ended with two taxa, having an inferred speciation rate of 0.
NULL