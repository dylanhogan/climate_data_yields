suppressMessages(here::i_am('R_analysis/figures/response_functions.R'))
CODE = here::here()
if (is.null(init)) source(glue::glue("{CODE}/R_analysis/init.R"))
local_load('analysis/regress.R')
local_load('analysis/clean_data.R')

generate_plot_df = function(
        df, source, crop,
        GDD_lower, KDD_lower,
        temp_lower, temp_upper,
        centering = 10, regression_args=list()) {
    
    # Initialize temperatures and output dataframe
    temp = data.frame(index = seq(1, 161*2, 1)) |> 
        mutate(temperature = (index - 1) / 8)
    output_df = data.frame()
    
    
    # Piecewise linear specification
    m = do.call(run_yield_regression, append(list(df=df, spec='piecewise'), regression_args))
    
    # Construct response function by evaluating at each temperature
    piecewise_hypothesis = function(t, m) {
        if (t <= GDD_lower) 
            h = hypotheses(m, hypothesis="0=0")
        else if ((t > GDD_lower) & (t <= KDD_lower))
            h = hypotheses(m, hypothesis=glue("GDD*{t-GDD_lower}=0"))
        else if ((t > KDD_lower) & (t <= temp_upper))
            h = hypotheses(m, hypothesis=glue("GDD*{KDD_lower - GDD_lower} + KDD*{t-KDD_lower}=0"))
        else if (t > temp_upper)
            h = hypotheses(m, hypothesis=glue("GDD*{KDD_lower - GDD_lower} + KDD*{temp_upper-KDD_lower}=0"))
        return(h)
    }
    
    # Generate dataframe
    pw_df = lapply(temp$temperature, piecewise_hypothesis, m=m) |> 
        bind_rows() |> 
        select(-term) |> 
        bind_cols(temp) |>
        select(temperature, b=estimate, lb=conf.low, ub=conf.high) |>
        mutate(
            lb = ifelse(temperature <= GDD_lower, 0, lb),
            ub = ifelse(temperature <= GDD_lower, 0, ub),
            spec = 'piecewise',
            clim = source,
            crop = crop
        )
    
    # Center response function
    center = pw_df |> filter(temperature == centering) |> unlist()
    pw_df = pw_df |>
        mutate(
            b = b - as.numeric(center['b']),
            lb = lb - as.numeric(center['b']),
            ub = ub - as.numeric(center['b'])
        )
    
    # Store dataframe
    output_df = bind_rows(output_df, pw_df)
    
    
    # Polynomial specification
    m = do.call(run_yield_regression, append(list(df=df, spec='poly'), regression_args))
    
    # Initialize Cheby terms
    poly = temp |>
        mutate(
            cheby = ifelse(
                temperature <=36,
                (temperature +3) / 39 * 2 - 1,
                NA
            ),
            cheby1 = cheby,
            cheby2 = 2*cheby^2 - 1
        )
    for (p in 3:8)
        poly = poly |> mutate(
            !! glue("cheby{p}") := 2 * cheby * !!sym(glue("cheby{p-1}")) - !!sym(glue("cheby{p-2}"))
        )
    poly = poly |> select(-cheby)
    
    C = poly |> 
        select(contains('cheby')) |> 
        as.matrix()
    
    # Evaluate polynomial at each temp
    cheby_hypothesis = function(i, m) {
        mult = \(x,y) paste(as.character(x), as.character(y), sep=' * ')
        hyp = paste(
            mapply(
                mult,
                x=names(C[i,]),
                y=as.character(C[i,]),
                SIMPLIFY=FALSE),
            collapse = " + "
        )
        return(hypotheses(m, hypothesis=glue("{hyp} = 0")))
    }
    
    # Construct output dataframe
    poly_df = lapply(temp$index, cheby_hypothesis, m=m) |> 
            bind_rows() |> 
            select(-term) |> 
            bind_cols(temp) 
    poly_df = poly_df |> 
            tidyr::replace_na(as.list(filter(poly_df, temperature==temp_upper))) |>
            select(temperature, b=estimate, lb=conf.low, ub=conf.high) |>
            mutate(spec='poly', clim=source, crop=crop)
    
    # Center response function
    center = poly_df |> filter(temperature == centering) |> unlist()
    poly_df = poly_df |>
        mutate(
            b = b - as.numeric(center['b']),
            lb = lb - as.numeric(center['b']),
            ub = ub - as.numeric(center['b'])
        )
    
    output_df = bind_rows(output_df, poly_df)
    
    
    # 3-degree bins
    m = do.call(run_yield_regression, append(list(df=df, spec='bins'), regression_args))
    
    # Evaluate at each bin
    bins = data.frame()
    j = \(x, t=0) gsub('-', 'Minus', x-t)
    epsilon = 0.001
    for (i in seq(0, temp_upper-3, 3)) {
        binv = glue::glue('bin_{j(i)}_{j(i+3)}')
        h = hypotheses(m, hypothesis=glue::glue('{binv}=0')) |>
            select(-term)
        bins = bind_rows(bins, h |> mutate(temperature = i + epsilon))
        bins = bind_rows(bins, h |> mutate(temperature = i+3 - epsilon))
        if ((centering >= i) & (centering < i+3)) cbin = i
    }
    
    h = hypotheses(m, hypothesis=glue::glue('bin_{j(temp_upper)}_Inf=0')) |>
        select(-term)
    bins = bind_rows(bins, h |> mutate(temperature = temp_upper + epsilon))
    bins = bind_rows(bins, h |> mutate(temperature = temp_upper+3 - epsilon))
    bins = bind_rows(bins, h |> mutate(temperature = 40))
    bins = bins |>
        select(temperature, b=estimate, lb=conf.low, ub=conf.high) |>
        mutate(spec = 'bins', clim = source, crop = crop)
    
    
    # Center bins
    center = bins |> filter(temperature == cbin + epsilon) |> unlist()
    bins = bins |> mutate(
        b = b - as.numeric(center['b']),
        lb = lb - as.numeric(center['b']),
        ub = ub - as.numeric(center['b']),
    )
    
    output_df = bind_rows(output_df, bins)
    
    return(output_df)
    
    
}

plot_response_function = function(spec, crop, plot_df, bounds = c(5, 40)) {
    
    temp_df = plot_df |> 
        filter(crop == !!crop, spec == !!spec, temperature > bounds[1], temperature <= bounds[2])

    temp_df$clim = factor(
        temp_df$clim,
        levels = c('prism', 'era5', 'gmfd'),
        labels = c('PRISM', 'ERA5-Land', 'GMFD')
    )
    temp_df$crop = factor(
        temp_df$crop,
        levels = c('corn', 'soy'),
        labels = c('Corn', 'Soybean')
    )
    temp_df$spec = factor(
        temp_df$spec,
        levels = c('piecewise', 'poly', 'bins', 'tavg'),
        labels = c('Piecewise linear', 'Polynomial', '3-degree bins', 'Average temp.')
    )
    spec_name = temp_df[1,'spec']
    crop_name = temp_df[1, 'crop']

    plot = ggplot(data=temp_df) +
        geom_ribbon(aes(ymin=lb, ymax=ub, x=temperature, fill=clim), alpha=.30) +
        geom_line(aes(y=b, x=temperature, color=clim), linewidth=.75) +
        scale_fill_manual(values=c('#EE6677', '#4477AA', '#228833')) +
        scale_color_manual(values=c('#EE6677', '#4477AA', '#228833')) +
        labs(
            color='Climate data set', fill='Climate data set',
            y=glue::glue('Log {crop_name} Yields (Bushel/Acre)'), x='Temperature (°C)',
            title=spec_name
        ) +
        theme_classic() + 
        scale_y_continuous(breaks=seq(-0.06, 0.02, 0.02), limits=c(-0.07, 0.02)) +
        scale_x_continuous(expand=c(0,0)) +
        theme(
            legend.position=c(.35, .35),
            legend.box.background = element_rect(colour = "black"),
            legend.text=element_text(size=14),
            text = element_text(size = 14),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line( linewidth=.1, color="gray"))

    if (spec != 'piecewise')
        plot = plot + theme(axis.title.y = element_blank())

    if (crop == 'corn') {
        plot = plot + theme(axis.title.x = element_blank())
    } else {
        plot = plot + theme(plot.title = element_blank())
    }
    return(plot)
}

plot_response_function_matrix = function(plot_df) {
    
    plots = mapply(
        plot_response_function,
        spec = c("piecewise", "poly", "bins"),
        crop = c(rep(c("corn"), 3), rep(c("soy"), 3)),
        MoreArgs = list(plot_df = plot_df),
        SIMPLIFY=FALSE
    )
    plot = ggarrange(
        plotlist=plots, ncol=3, nrow=2,
        legend='bottom', common.legend=TRUE
    )
    return(plot)
    
}

generate_response_function_plot_df = function(
    year_min, year_max,
    growing_season_min, growing_season_max,
    GDD_lower, KDD_lower_list,
    temp_lower, temp_upper,
    centering = 10, regression_args=list(),
    sources = c("prism", "era5", "gmfd"),
    crops = c("corn", "soy")
){
    
    df_list = list()
    for (source in sources) {
        for (crop in crops) {
            
            KDD_lower = KDD_lower_list[[crop]][[source]]
            
            # Load standard panel data set
            df = load_panel(
                source, crop,
                year_min, year_max,
                growing_season_min, growing_season_max,
                GDD_lower, KDD_lower,
                temp_lower, temp_upper
            )
            
            # Run regressions and create DF of evaluated responses
            plot_df = generate_plot_df(
                df, source, crop, 
                GDD_lower, KDD_lower,
                temp_lower, temp_upper,
                centering, regression_args
            )
            df_list = append(df_list, list(plot_df))
        }
    }
    
    return(bind_rows(df_list))
}

temperature_sensitivity_table = function(
    year_min, year_max,
    growing_season_min, growing_season_max,
    GDD_lower, KDD_lower_list,
    temp_lower, temp_upper,
    centering = 10, regression_args=list(),
    sources = c("prism", "era5", "gmfd"),
    crops = c("corn", "soy")
) {

    impact_df = list()
    for (source in sources) {
        for (crop in crops) {

            KDD_lower = KDD_lower_list[[crop]][[source]]

            # Load standard panel data set
            df = load_panel(
                source, crop,
                year_min, year_max,
                growing_season_min, growing_season_max,
                GDD_lower, KDD_lower,
                temp_lower, temp_upper
            )        
            
            # Run regression and calculate impacts with uncertainty
            m = do.call(run_yield_regression, append(list(df=df, spec='bins'), regression_args))

            ref_bin = 'bin_9_12'
            hot_bin = 'bin_36_Inf'
            max_bin = names(which.max(coeftable(m)[,'Estimate']))

            hot_impact = hypotheses(m, hypothesis=glue('({hot_bin} - {ref_bin})*100 = 0'))
            max_impact = hypotheses(m, hypothesis=glue('({max_bin} - {ref_bin})*100 = 0'))

            temp_impact = bind_rows(hot_impact, max_impact) |> 
                select(b=estimate, se=std.error) |>
                mutate(
                    crop = !!crop,
                    clim = !!source,
                    temp = c('hot', 'warm')
                )

            impact_df = bind_rows(impact_df, temp_impact)

        }
    }
    
    # Generate table
    anon = function(x) {paste(glue('({sprintf("%.2f", x)})'))}
    tab = data.table::melt(data.table::as.data.table(impact_df), id.vars=c('clim','crop','temp')) |>
        mutate(value = round(value, 2)) |>
        spread(variable, value) |>
        mutate(
            se = anon(se),
            b = sprintf("%.2f", b)
        ) |>
        gather(variable, value, b , se) |>
        spread(clim, value) |>
        select(crop, temp, variable, prism, era5, gmfd)

    return(xtable(tab))
}