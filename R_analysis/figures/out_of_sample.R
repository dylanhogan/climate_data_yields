suppressMessages(here::i_am('R_analysis/figures/out_of_sample.R'))
CODE = here::here()
if (is.null(init)) source(glue::glue("{CODE}/R_analysis/init.R"))
local_load('analysis/regress.R')
local_load('analysis/clean_data.R')

generate_out_of_sample_plot_df = function(
    year_min, year_max,
    growing_season_min, growing_season_max,
    GDD_lower, KDD_lower_list,
    temp_lower, temp_upper,
    regression_args=list(),
    num_iterations=1000, num_cores=30,
    sources = c("prism", "era5", "gmfd"),
    crops = c("corn", "soy"),
    transform = c('piecewise', 'bins', 'poly', 'tavg'),
    seed=123
) {
    
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
                temp_lower, temp_upper,
                transform
            )
            
            # Run out of sample prediction exercise
            temp = run_out_of_sample_prediction(
                df,
                spec_list = c('baseline', transform),
                num_iterations=num_iterations,
                num_cores=num_cores,
                regression_args=regression_args,
                seed=seed
            ) |>
            mutate(clim = !!source, crop = !!crop)

            df_list = append(df_list, list(temp))
        }
    }
    
    oos = bind_rows(df_list)
    
    return(oos)
    
}

plot_rms = function(
    oos, 
    transform = c('piecewise', 'bins', 'poly', 'tavg'),
    colors = c('#EE6677', '#4477AA', '#228833')
) {
    
    bar = plot_rms_bar(oos, transform, colors)
    box = plot_rms_box(oos, transform, colors)
    
    leg = ggpubr::get_legend(bar)
    return (
        cowplot::plot_grid(box, leg, ncol=1, rel_heights = c(15, 1))
    )
}

plot_rms_bar = function(
    oos, 
    transform = c('piecewise', 'bins', 'poly', 'tavg'),
    colors = c('#EE6677', '#4477AA', '#228833')
 ) {
    
    oos = oos |>
        select(-iteration) |>
        group_by(clim, crop) |>
        summarize_all(mean) 
    for (spec in transform) {
        oos = oos|>
            mutate(!! spec := (1 - !!sym(spec)/baseline)*100)
    }
    oos = oos |> select(-baseline)
    
    plot_df = melt(as.data.table(oos), id.vars=c('clim','crop')) |> 
        rename(spec=variable, rmse=value)
    
    plot_df$clim = factor(
        plot_df$clim,
        levels = c('prism', 'prism_resample', 'era5', 'era5_sine', 'gmfd'),
        labels = c(
            'PRISM', 'PRISM\nResampled to ERA5 grid', 
            'ERA5-Land', 'ERA5-Land\nSine-interpolated', 
            'GMFD'
        )
    )
    plot_df$crop = factor(
        plot_df$crop,
        levels = c('corn', 'soy'),
        labels = c('Corn', 'Soybeans')
    )
    plot_df$spec = factor(
        plot_df$spec,
        levels = c('piecewise', 'poly', 'bins', 'tavg'),
        labels = c('Piecewise linear', 'Polynomial', '3-degree bins', 'Average temp.')
    )
    
    plot = ggplot(plot_df, aes(x=clim, y = rmse, fill=clim, alpha=spec)) + 
        geom_col(position="dodge", color='black') +
        scale_fill_manual(values = colors, guide="none") +
        scale_alpha_manual(values = c(.9,0.7,0.5,0.3),name="") +
        scale_y_continuous(expand=c(0,0), breaks=c(0,5,10,15)) +
        facet_wrap(~crop) +
        theme_classic() +
        labs( y='Percent reduction in RMS', alpha='') +
        theme(
            legend.position='bottom',
            legend.box.background = element_rect(colour = "black"),
            axis.text=element_text(size=14),
            axis.title=element_text(size=14),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line( linewidth=.1, color="grey"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=14),
            strip.background = element_blank(),
            strip.text = element_text(size=14),
            legend.text=element_text(size=14)
        ) 
    
    return(plot)
}

plot_rms_box = function(
        oos, 
        transform = c('piecewise', 'bins', 'poly', 'tavg'),
        colors = c('#EE6677', '#4477AA', '#228833')) {
    
    plot_df = oos
    for (spec in transform) {
        plot_df = plot_df |>
            mutate(!! spec := (1 - !!sym(spec)/baseline)*100)
    }
    plot_df = plot_df |> select(-baseline)
    plot_df = pivot_longer(
        plot_df, 
        cols=c(piecewise, bins, poly, tavg), 
        names_to='spec',
        values_to='rmse'
    )
    
    means = plot_df |> 
        group_by(crop, clim, spec) |> 
        summarize(mean_rmse = mean(rmse, na.rm=TRUE), .groups='drop') 
    plot_df = plot_df |> left_join(means, by=c('clim', 'crop', 'spec'))
    
    plot_df$clim = factor(
        plot_df$clim,
        levels = c('prism', 'prism_resample', 'era5', 'era5_sine', 'gmfd'),
        labels = c(
            'PRISM', 'PRISM\nResampled to ERA5 grid', 
            'ERA5-Land', 'ERA5-Land\nSine-interpolated', 
            'GMFD'
        )
    )
    plot_df$crop = factor(
        plot_df$crop,
        levels = c('corn', 'soy'),
        labels = c('Corn', 'Soybeans')
    )
    plot_df$spec = factor(
        plot_df$spec,
        levels = c('piecewise', 'poly', 'bins', 'tavg'),
        labels = c('Piecewise linear', 'Polynomial', '3-degree bins', 'Average temp.')
    )
    
    dodge <- position_dodge(width=0.75) 
    plot = ggplot(plot_df, aes(x=clim, y = rmse, fill=clim, alpha=spec, coef = 0)) + 
        geom_boxplot(outlier.shape=NA, coef=0) +
        geom_point(aes(y=mean_rmse), position=dodge, size=1.5) +
        scale_fill_manual(values = colors, guide='none') +
        scale_alpha_manual(values = c(.9,0.7,0.5,0.3), guide='none') +
        scale_y_continuous(expand=c(0,0), breaks=c(0,5,10,15, 20)) +
        facet_wrap(~crop) +
        theme_classic() + 
        coord_cartesian(ylim = c(0,19)) +
        labs( y='Percent reduction in RMS') +
        theme(
            legend.position='bottom',
            legend.box.background = element_rect(colour = "black"),
            axis.text=element_text(size=14),
            axis.title=element_text(size=14),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line( linewidth=.1, color="grey"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=14),
            strip.background = element_blank(),
            strip.text = element_text(size=14),
            legend.text=element_text(size=14)
        ) 
    
    return(plot)
}