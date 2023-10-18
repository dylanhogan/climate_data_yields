suppressMessages(here::i_am('R_analysis/analysis/clean_data.R'))
CODE = here::here()
if (is.null(init)) source(glue("{CODE}/R_analysis/init.R"))

load_climate_data = function(
    source, crop,
    year_min=1950, year_max = 2020,
    growing_season_min = 3, growing_season_max = 9
) {

    pth = glue("{DATA}/{source}/weather_{crop}.feather")
    
    if (source %in% c('era5', 'gmfd', 'prism')) {
        # Open data set, filter to years/months
        weather = read_feather(pth) 
        weather = weather |>
            select(fips, year, month, contains("dday"), contains("time"), tavg, prec) |>
            filter(year >= year_min, year <= year_max) |>
            filter(month >= growing_season_min, month <= growing_season_max) |>
            select(-month)
        
        # Aggregate over growing season
        weather = mapply(summarise_at, 
                .vars = lst(setdiff(names(weather), c("fips", "year", "tavg")), "tavg"), 
                .funs = lst(sum, mean), 
                MoreArgs = list(.tbl = weather |> group_by(fips, year))) |>
            reduce(merge, by = c("fips", "year"))
        
        weather = weather |>
            mutate(prec2 = prec^2, tavg2 = tavg^2)
        
    } else {
        # Open data set, filter to years/months
        weather = read_stata(pth) 
        weather = weather |>
            select(fips, year, month, contains("dday"), prec) |>
            filter(year >= year_min, year <= year_max) |>
            filter(month >= growing_season_min, month <= growing_season_max) |>
            select(-month)
        
        # Aggregate over growing season
        weather = weather |> 
            group_by(fips, year) |>
            summarize_all(mean) |>
            ungroup()
        
        weather = weather |>
            mutate(prec2 = prec^2)
        
    }
    
    # Standardize precip units
    if (source == 'era5') weather = mutate(weather, prec = prec * 1000)
    if (source == 'gmfd') weather = mutate(weather, prec = prec * 100000)
    
    return(weather)
    
}

load_yield_data = function(crop, year_min=1950, year_max = 2020) {
    cn = crop
    crops = read_feather(glue("{DATA}/yields/yieldData.feather")) |>
        as.data.frame() |>
        filter(crop == cn, year >= year_min, year <= year_max)
    return(crops)
}


cheby_poly8 = function(df, temp_lower = -3, temp_upper = 36) {
    
    # Generate Cheby matrix
    C = matrix(nrow=temp_upper-temp_lower, ncol=8)
    for (t in seq(temp_lower, temp_upper-1)) {
        ind = t-temp_lower + 1
        C[ind, 1] = (t + 0.5 - temp_lower) / (temp_upper - temp_lower)*2 - 1
        C[ind, 2] = 2 * C[ind,1] * C[ind,1] - 1
        for (c in 3:8) C[ind, c] = 2 * C[ind,1] * C[ind,c-1] - C[ind,c-2]
    }
    
    # Generate Cheby vars
    for (p in 1:8) {
        varn = glue('cheby{p}')
        df[varn] = df[glue('time{temp_upper}C')] * C[temp_upper-1-temp_lower+1,p]
        for (t in seq(temp_lower, temp_upper-2)) {
            lb = gsub('-', 'Minus', t)
            ub = gsub('-', 'Minus', t+1)
            ind = t-temp_lower + 1
            df[varn] = df[varn] + (df[glue('time{lb}C')] - df[glue('time{ub}C')])*C[ind,p]
        }
    }
    
    return(df)
    
}

time_bin3 = function(df, temp_lower = -3, temp_upper = 36) {
    
    # Generate bins
    for (t in seq(temp_lower, temp_upper-3, 3)) {
        lb = gsub('-', 'Minus', t)
        ub = gsub('-', 'Minus', t+3)
        varn = glue('bin_{lb}_{ub}')
        df[varn] = df[glue('time{lb}C')] - df[glue('time{ub}C')]
    }
    varn = glue('bin_{temp_upper}_Inf')
    df[varn] = df[glue('time{temp_upper}C')]
    
    return(df)
    
}

piecewise = function(
        df, GDD_lower, KDD_lower,
        temp_lower = -3, temp_upper = 36
    ) {

    # Define strings for piecewise linear spec kinks
    d = \(x, t=0) sym(gsub('-', 'Minus', glue('dday{x-t}C')))
    
    # Generate GDD, KDD vars
    df = df |>
        mutate(
            KDD = !!d(KDD_lower) - !!d(temp_upper),
            GDD = !!d(GDD_lower) - !!d(KDD_lower)
        )
    
    return(df)
}

transform_temperature = function(
    df, GDD_lower=10, KDD_lower=NULL,
    temp_lower = -3, temp_upper = 36,
    transform = c('piecewise', 'bins', 'poly', 'tavg')
) {
    
    if ('piecewise' %in% transform)
        df = piecewise(df, GDD_lower, KDD_lower, temp_lower, temp_upper)
    if ('bins' %in% transform)
        df = time_bin3(df, temp_lower, temp_upper)
    if ('poly' %in% transform)
        df = cheby_poly8(df, temp_lower, temp_upper)
    
    return(df)
}

clean_panel = function(
    panel, GDD_lower, KDD_lower, 
    temp_lower, temp_upper,
    transform = c('piecewise', 'bins', 'poly', 'tavg')
) {

    # Generate variables for analysis
    panel = panel |>
        mutate(
            log_yield = log(yield),
            state = floor(fips/1000),
            t = year - min(year),
            t2 = t^2,
        ) 
    
    panel = transform_temperature(
            panel, GDD_lower, KDD_lower,
            temp_lower, temp_upper, transform
        )
    return(panel)
}


load_panel = function(
    source, crop,
    year_min=1950, year_max = 2020,
    growing_season_min = 3, growing_season_max = 9,
    GDD_lower=10, KDD_lower=30, 
    temp_lower=-3, temp_upper=36,
    transform = c('piecewise', 'bins', 'poly', 'tavg')
) {
    
    # Load weather
    weather = load_climate_data(
        source, crop, year_min, year_max,
        growing_season_min, growing_season_max
    )
    
    # Load crops
    crops = load_yield_data(crop, year_min, year_max)
    
    # Merge panel
    panel = inner_join(crops, weather, by=c('fips', 'year'))
    
    # Clean panel
    panel = clean_panel(panel,
        GDD_lower, KDD_lower,
        temp_lower, temp_upper, transform)

    return(panel)
    
}
load_panel = memoise(load_panel, cache=cm)

load_main_analysis_panel = function(source, crop) {
    
    # Load arguments for main analysis
    source(glue("{CODE}/R_analysis/analysis/main_analysis_args.R"))
    
    df = load_panel(
        source, crop,
        year_min, year_max,
        growing_season_min, growing_season_max,
        GDD_lower, KDD_lower,
        temp_lower, temp_upper
    )
    
    return(df)

}
    
    
    