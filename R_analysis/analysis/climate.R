suppressMessages(here::i_am('R_analysis/analysis/climate.R'))
CODE = here::here()
if (is.null(init)) source(glue("{CODE}/R_analysis/init.R"))
local_load('analysis/regress.R')

piecewise_climate_shift = function(df, t, b) {
    yhat = glue::glue('yhat{t}')
    prod = glue::glue('prod{t}')
    KDD_v = sym(glue::glue('KDD_clim{t}'))
    GDD_v = sym(glue::glue('GDD_clim{t}'))
    df = df |> mutate(
        !! yhat := yhat0 + 
            b['KDD'] * !!KDD_v +
            b['GDD'] * !!GDD_v,
        !! prod := exp(!!sym(yhat))*areaHarv
    )
    return(df)
}

bins_climate_shift = function(df, t, b) {
    yhat = glue::glue('yhat{t}')
    prod = glue::glue('prod{t}')
    
    df = df |> mutate(!! yhat := yhat0)
    bins = names(df)[
        grepl('bin', names(df)) &
        grepl(glue::glue('clim{t}'), names(df)) &
        ! grepl('Inf', names(df))
    ]
    
    for (bin in bins) {
        varn = gsub(glue::glue('_clim{t}'),'',bin)
        varn_shift = sym(bin)
        df = df |> mutate(!! yhat := !!sym(yhat) + b[varn] * !!varn_shift)
    }
    varn = gsub(
        glue::glue('_clim{t}'), '',
        names(df)[grep(glue::glue('Inf_clim{t}'), names(df))]
    )
    varn_shift = sym(names(df)[grep(glue::glue('Inf_clim{t}'), names(df))])
    df = df |> mutate(
        !! yhat := !!sym(yhat) + b[varn] * !!varn_shift,
        !! prod := exp(!!sym(yhat))*areaHarv
    )
    return(df)
}


poly_climate_shift = function(df, t, b) {
    yhat = glue::glue('yhat{t}')
    prod = glue::glue('prod{t}')
    
    df = df |> mutate(!! yhat := yhat0)
    for (p in 1:8) {
        varn = glue::glue('cheby{p}')
        varn_shift = sym(glue::glue('{varn}_clim{t}'))
        df = df |> mutate(!! yhat := !!sym(yhat) + b[varn] * !!varn_shift)
    }
    df = df |> mutate(!! prod := exp(!!sym(yhat))*areaHarv)
    return(df)
}


tavg_climate_shift = function(df, t, b) {
    yhat = glue::glue('yhat{t}')
    prod = glue::glue('prod{t}')
    tavg_v = sym(glue::glue('tavg_clim{t}'))
    tavg2_v = sym(glue::glue('tavg2_clim{t}'))
    df = df |> mutate(
        !! yhat := yhat0 + 
            b['tavg'] * !!tavg_v +
            b['tavg2'] * !!tavg2_v,
        !! prod := exp(!!sym(yhat))*areaHarv
    )
    return(df)
}


generate_shifted_climate_variables = function(
    df, GDD_lower, KDD_lower, temp_upper, temp_lower
){
    
    # Generate Cheby matrix
    C = matrix(nrow=temp_upper-temp_lower, ncol=8)
    for (t in seq(temp_lower, temp_upper-1)) {
        ind = t-temp_lower + 1
        C[ind, 1] = (t + 0.5 - temp_lower) / (temp_upper - temp_lower)*2 - 1
        C[ind, 2] = 2 * C[ind,1] * C[ind,1] - 1
        for (c in 3:8) C[ind, c] = 2 * C[ind,1] * C[ind,c-1] - C[ind,c-2]
    }
    
    d = \(x, t=0) sym(gsub('-', 'Minus', glue::glue('dday{x-t}C')))
    j = \(x, t=0) gsub('-', 'Minus', max(x-t, -4))
    for (t in 1:4) {
        
        # Piecewise linear
        df = df |> mutate(
            !! glue::glue('KDD_clim{t}') := !!d(KDD_lower,t) - !!d(temp_upper,t) - !!d(KDD_lower) + !!d(temp_upper),
            !! glue::glue('GDD_clim{t}') := !!d(GDD_lower,t) - !!d(KDD_lower,t) - !!d(GDD_lower) + !!d(KDD_lower)
        )
        
        # 3-degree bins
        # Generate shifted bins
        for (i in seq(temp_lower, temp_upper-3, 3)) {
            varn = glue::glue('bin_{j(i)}_{j(i+3)}')
            varn_shift = glue::glue('{varn}_clim{t}')
            # print(glue::glue("time{j(i,t)}C"))
            df[varn_shift] = (df[glue::glue('time{j(i,t)}C')] - df[glue::glue('time{j(i+3,t)}C')]) - df[varn]
        }        
        varn = glue::glue('bin_{j(temp_upper)}_Inf')
        varn_shift = glue::glue('{varn}_clim{t}')
        df[varn_shift] = (df[glue::glue('time{j(temp_upper,t)}C')]) - df[varn]
        
        # Generate Cheby vars
        for (p in 1:8) {
            varn = glue::glue('cheby{p}_clim{t}')
            df[varn] = df[glue::glue('time{j(temp_upper,t)}C')] * C[temp_upper-1-temp_lower+1,p]
            for (i in seq(temp_lower, temp_upper-2)) {
                # lb = gsub('-', 'Minus', t)
                # ub = gsub('-', 'Minus', t+1)
                ind = i - temp_lower + 1
                df[varn] = df[varn] + (df[glue::glue('time{j(i,t)}C')] - df[glue::glue('time{j(i+1,t)}C')])*C[ind,p]
            }
            df[varn] = df[varn] - df[glue::glue('cheby{p}')]
        }
        
        # Average temp.
        varn = glue::glue('tavg_clim{t}')
        varn2 = glue::glue('tavg2_clim{t}')
        tv = t
        df = df |> mutate(
            !! varn := tv,
            !! varn2 := (tavg + tv)^2 - tavg^2
        )
        
    }
    return(df)
}

climate_impact = function(df, b, shift_func, t_list=c(1,2,3,4), aggregate=TRUE) {
    
    for (t in t_list) df = df |> shift_func(t, b)

    pdf = df |>
        select(fips, year, contains('prod')) 
    
    if (aggregate) {
        pdf = pdf |>
            select(-fips) |>
            group_by(year) |>
            summarize_all(sum) |>
            select(-year) |>
            summarize_all(mean)
    } else {
        pdf = pdf |>
            select(-year) |>
            group_by(fips) |>
            summarize_all(mean)
    }

    for (t in t_list) {
        prod = sym(glue::glue('prod{t}'))
        impact = glue::glue('impact{t}')
        pdf = pdf |> 
            mutate(!!impact := (!!prod / prod0) - 1) 
    }
    pdf = pdf |> select(-contains('prod'))
    return(pdf)
    
}


run_climate_projection = function(
    df,
    spec,
    GDD_lower, KDD_lower, 
    temp_lower, temp_upper,
    aggregate=TRUE, t_list=c(1,2,3,4),
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
    
    m = do.call(run_yield_regression, append(list(df=df, spec=spec), regression_args))
    
    # Store residuals
    df$residuals = resid(m)
    df$yhat0 = predict(m, df)
    df = df |> 
        dplyr::filter(year >=1960, year < 1989) |>
        mutate(prod0 = exp(yhat0)*areaHarv) |>
        generate_shifted_climate_variables(GDD_lower, KDD_lower, temp_upper, temp_lower) |> 
        select(crop, fips, year, yhat0, prod0, residuals, areaHarv, contains('clim'))
    
    # Point estimate
    vars = get_temp_vars(spec)
    b = coef(m)[vars]
    v = vcov(m)[vars, vars]
    
    shift_str = glue::glue('{spec}_climate_shift')
    shift_func = get(shift_str)
    
    point = climate_impact(df, b, shift_func, t_list, aggregate)
    if (aggregate) point$iteration = 0
    
    # Perform the process in parallel and store the impact from each iteration
    if (num_iterations > 0) {
        
        # Create matrix of draws from uncertainty
        draws_mat = MASS::mvrnorm(num_iterations, mu=b, Sigma=v)
        
        clim_results = foreach(
            i = 1:num_iterations, 
            .export = c(shift_str, 'climate_impact'),
            .packages = c('dplyr', 'fixest')
        ) %dorng% {

            # Take draw
            b = draws_mat[i,]

            # Calculate impact
            impact = climate_impact(
                df, b, shift_func, 
                t_list, aggregate
            )

            # Return impact
            return(impact)
        }

        # Combine the results into one dataframe
        if (aggregate) {
            out = bind_cols(
                data.frame(iteration=1:num_iterations), 
                bind_rows(clim_results) 
            )
        } else {
            out = bind_rows(clim_results)
        }
    } else {
        out = data.frame()
    }
    
    # Including point estimate
    out = bind_rows(point, out) |>
        mutate(spec=!!spec)
    
    return(out)
}



