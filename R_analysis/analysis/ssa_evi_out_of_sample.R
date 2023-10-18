suppressMessages(here::i_am('R_analysis/analysis/ssa_evi_out_of_sample.R'))
CODE = here::here()
source(glue::glue("{CODE}/R_analysis/init.R"))

load_evi_panel = function(clim) {
    
    if (clim == "era5") {
        weather = glue("{DATA}/era5/era5_ssa_evi.feather")
        wdf = read_feather(weather) |>
            mutate(
                GDD = corn_dday10 - corn_dday27,
                KDD = corn_dday27 - corn_dday36
            )
    } else if (clim == "gmfd") {
        weather = glue("{DATA}/gmfd/gmfd_ssa_evi.feather")
        wdf = read_feather(weather) |>
            mutate(
                GDD = gmfd_corn_dday10 - gmfd_corn_dday28,
                KDD = gmfd_corn_dday28 - gmfd_corn_dday36
            )
    } else if (clim == "cru") {
        weather = glue("{DATA}/cru/cru_ssa_evi.feather")
        wdf = read_feather(weather) |>
            mutate(
                GDD = cru_corn_dday10 - cru_corn_dday25,
                KDD = cru_corn_dday25 - cru_corn_dday36
            ) |>
            rename(pr=pre)
    }
    
    wdf = wdf |> filter(year <= 2010)

    
    crop_stations = read_feather(glue('{DATA}/evi/ssa_station_and_crop.feather'))
    crop_stations$station = 1
    wdf = right_join(crop_stations, wdf, by=c('longitude', 'latitude')) |>
        mutate(station = ifelse(is.na(station), 0, 1))
    
    evi = read_feather(glue("{DATA}/evi/ssa_evi.feather"))
    df = evi |> 
        inner_join(wdf, by=c('latitude', 'longitude', 'year', 'month'))

    gs = read_csv(glue('{DATA}/evi/growing_seasons.csv'), show_col_types = FALSE)
    df = df |>  
        rename(region=ISO) |>
        left_join(gs, by='region')
    df = df |>
        mutate(
            skip_year = (month_start > month_end),
            is_gs = ifelse(
                month_start > month_end,
                (month <= month_end) | (month >= month_start),
                (month >= month_start) & (month <= month_end)
            )
        ) |>
        filter(is_gs) |>
        mutate(
            year = ifelse(
                skip_year,
                ifelse(
                    month >= month_start,
                    year+1,
                    year
                ),
                year
            )
        ) |>
        select(-month_start, -month_end, -skip_year, -is_gs) |>
        arrange(region, year, month)
    
    df = df |> 
        arrange(latitude, longitude, year, month) |> 
        group_by(latitude, longitude) |> 
        mutate(id = cur_group_id())
    
    estim = df |> 
        group_by(region, id, year) |>
        summarize(
            GDD = sum(GDD),
            KDD = sum(KDD),
            pr = sum(pr),
            total_evi = sum(evi/10000),
            max_evi = max(evi/10000),
            station = mean(station)
        ) |>
        ungroup() |>
        mutate(
            log_total_evi = log(total_evi),
            log_max_evi = log(max_evi),
            pr_2 = pr^2,
            t = year - min(year),
            t_2 = t^2
        ) |> 
        group_by(region) |> 
        mutate(region_id = cur_group_id()) |>
        filter(id != 100000, log_total_evi > 0)

    estim = estim |> 
        filter(year < 2010, year > 2000)
    
   return(estim)
    
}
load_evi_panel = memoise(load_evi_panel, cache=cm)

run_evi_out_of_sample = function(
        df, 
        station_filters=c('only_stations', 'only_not_stations'),
        seed=123, num_iterations=1000, num_cores=20) {

    # Create a parallel backend for parallel execution
    cl = makeCluster(num_cores)

    # Register the parallel backend
    registerDoParallel(cl)
    
    # Set the seed for reproducibility
    set.seed(seed, kind = "L'Ecuyer-CMRG")
                      
    sample_years = function(years) {
        # Split the data into train and test sets based on the year
        train_indices = sample(length(years), floor(0.7 * length(years)))
        train_years = years[train_indices]
        test_years = years[-train_indices]
        return(list(train=train_years, test=test_years))
    }                      
                      
    # Create a function to perform the regression, prediction, and RMSE calculation
    regression_prediction_rmse = function(train_data, test_data, baseline=FALSE) {
        library(fixest)
        
        # Perform the regression on the train set
        if (baseline)
             model_degree_day = as.formula("log_total_evi ~ 1 | id + region_id[t] + region_id[t_2]")
        else
            model_degree_day = as.formula("log_total_evi ~ GDD + KDD + pr + pr_2 | id + region_id[t] + region_id[t_2]")
        
        model = feols(model_degree_day, data=train_data, cluster=~region_id)

        # Make predictions on the test set
        predicted_values = predict(model, newdata = test_data)

        # Calculate RMSE for the out-of-sample predictions
        actual_values = test_data$log_total_evi 

        rmse = sqrt(mean((actual_values - predicted_values)^2))

        # Return the RMSE
        return(rmse)
    }
    some_stations = df |>
        select(id, region, station) |>
        distinct() |>
        group_by(region) |> 
        summarize(station = sum(station))
    # Perform the process in parallel and store the RMSE for each iteration
    rmse_results = foreach(i = 1:num_iterations, .packages = "dplyr") %dorng% {
        
        rmse = data.frame()
        for (crop_station_filter in station_filters) {
            
            
            if (crop_station_filter != 'all') {
                if (crop_station_filter == 'only_stations') {
                    estim = filter(df, station==1)
                }
                else if (crop_station_filter == 'only_not_stations') {
                    no_stations = filter(df, station != 1)
                    # stations = unique(filter(df, station==1)$id)
                    estim = data.frame()
                    for (region in some_stations$region) {
                        stations = filter(some_stations, region == !!region)$station
                        no_station_ids = unique(filter(no_stations, region==!!region)$id)
                        stdraws = sample(no_station_ids, min(stations,length(no_station_ids)))
                        estim = bind_rows(estim, filter(no_stations, id %in% stdraws))
                    }
                }
            } else {
                estim = df
                stations = unique(filter(df, station==1)$id)
                stdraws = sample(estim$id, length(stations))
                estim = estim |> filter(id %in% stdraws)
            }
            
            y = sample_years(unique(df$year))
            
            train_data = estim |> dplyr::filter(year %in% y[['train']])
            test_data = estim |> dplyr::filter(year %in% y[['test']])
            
            # Split the data, run regression, predict, and calculate RMSE
            rmse = bind_rows(rmse,
                data.frame(
                    iteration=i,
                    rmse = regression_prediction_rmse(train_data, test_data),
                    baseline = regression_prediction_rmse(train_data, test_data, baseline=TRUE),
                    crop_station = crop_station_filter,
                    n = nrow(estim)
                )
            )
        }

        # Return the RMSE
        return(rmse)
    }
    
    out = bind_rows(rmse_results)
    
    stopCluster(cl)
    
    return(out)
                             
}
