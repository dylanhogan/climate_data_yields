out_of_sample_prediction = function(df, spec, num_iterations = 1000, num_cores = 20, seed = 123) {
    
    # Set RNG seed
    set.seed(seed)

    regression_prediction_rmse = function(iter, data, spec) {

        # Split the data into train and test sets based on the year
        years = unique(data$year)
        train_indices = sample(length(years), floor(0.85 * length(years)))
        train_years = years[train_indices]
        test_years = years[-train_indices]

        train_data = data |> dplyr::filter(year %in% train_years)
        test_data = data |> dplyr::filter(year %in% test_years)

        # Perform the regression on the train set
        model = train_data |> run_yield_regression(spec)

        # Make predictions on the test set
        predicted_values = predict(model, newdata = test_data)
        # predicted_values = predicted_values[!is.na(predicted_values)]

        # Calculate RMSE for the out-of-sample predictions
        actual_values = test_data$log_yield  # Replace 'your_dependent_variable' with the actual column name
        # actual_values = actual_values[!is.na(predicted_values)]

        rmse = sqrt(mean((actual_values - predicted_values)^2, na.rm=T))

        # Return the RMSE
        return(rmse)
    }
    
    print(regression_prediction_rmse(df, spec))
    
    # Perform the process in parallel and store the RMSE for each iteration
    rmse_results = foreach(i = 1:num_iterations, .export='FUN') %dopar% {
        
        
        # Split the data, run regression, predict, and calculate RMSE
        rmse = regression_prediction_rmse(df, spec)

        # Return the RMSE
        return(rmse)
    }

    # Combine the results into one dataframe
    results_df = data.frame(iteration = 1:num_iterations, RMSE = unlist(rmse_results))
    results_df$clim = clim

    results_list[[clim]] = results_df
    
    return(results_df)
}