suppressMessages(here::i_am('R_analysis/analysis/regress.R'))
CODE = here::here()
if (is.null(init)) source(glue("{CODE}/R_analysis/init.R"))

get_temp_vars = function(spec) {
    
    # Define temperature variables for various specifications
    if (spec == 'piecewise') 
        f_T = c("GDD","KDD")
    else if (spec == "bins")
        f_T = c('bin_Minus3_0', 'bin_0_3', 'bin_3_6', 
            'bin_6_9', 'bin_9_12', 'bin_12_15', 
            'bin_15_18', 'bin_18_21', 'bin_21_24', 
            'bin_24_27', 'bin_27_30', 'bin_30_33', 
            'bin_33_36', 'bin_36_Inf')
    else if (spec == 'poly')
        f_T = c('cheby1','cheby2','cheby3','cheby4',
                'cheby5','cheby6','cheby7','cheby8')
    else if (spec == 'tavg') 
        f_T = c("tavg", "tavg2")
    else
        f_T = ''
    
    return(f_T)
}

run_yield_regression = function(
    df, spec,
    g_P = 'prec + prec2',
    controls = 'fips + state[t] + state[t2]',
    cluster = 'fips + state',
    residual_dir = NULL
) {
    
    f_T = paste(get_temp_vars(spec), collapse=" + ")
    
    # Setup formula for regression
    model_degree_day = as.formula(glue::glue("log_yield ~ {f_T} + {g_P} | {controls}"))
    
    if (spec == 'baseline')
        model_degree_day = as.formula(glue::glue("log_yield ~ 1 | {controls}"))

    # Run regression, cluster SEs
    mtemp = fixest::feols(model_degree_day, data=df)
    
    if (is.list(cluster))
        m = summary(mtemp, vcov = cluster)
    else
        m = summary(mtemp, cluster=cluster)
    
    return(m)
}


run_out_of_sample_prediction = function(
    df,
    spec_list = c('baseline', 'piecewise', 'bins', 'poly', 'tavg'), 
    num_iterations = 1000, 
    num_cores = 20, 
    seed = 123,
    regression_args=list()
) {
    # Set RNG seed
    set.seed(seed)
    
    # Create a parallel backend for parallel execution
    cl = makeCluster(num_cores)

    # Register the parallel backend
    registerDoParallel(cl)

    sample_years = function(years) {
        # Split the data into train and test sets based on the year
        train_indices = sample(length(years), floor(0.85 * length(years)))
        train_years = years[train_indices]
        test_years = years[-train_indices]
        return(list(train=train_years, test=test_years))
    }
    
    regression_prediction_rmse = function(spec, train_data, test_data, regression_args) {
        
        # Perform the regression on the train set
        model = do.call(run_yield_regression, append(list(df=train_data, spec=spec), regression_args))

        # Make predictions on the test set
        predicted_values = predict(model, newdata = test_data)

        # Calculate RMSE for the out-of-sample predictions
        actual_values = test_data$log_yield 
        rmse = sqrt(mean((actual_values - predicted_values)^2, na.rm=T))
        
        # Return the RMSE
        return(rmse)
    }
    
    # Perform the process in parallel and store the RMSE for each iteration
    rmse_results = foreach(
        i = 1:num_iterations, 
        .export = c("run_yield_regression", "get_temp_vars")
    ) %dorng% {
        
        y = sample_years(unique(df$year))

        train_data = df |> dplyr::filter(year %in% y[['train']])
        test_data = df |> dplyr::filter(year %in% y[['test']])
        
        rmse = mapply(
            regression_prediction_rmse,
            spec=spec_list,
            MoreArgs = list(
                train_data=train_data,
                test_data=test_data,
                regression_args=regression_args
            )
        )
        
        # Return the RMSE
        return(rmse)
    }

    # Combine the results into one dataframe
    out = bind_cols(
        data.frame(iteration=1:num_iterations), 
        bind_rows(rmse_results) 
    )
    
    return(out)
}

