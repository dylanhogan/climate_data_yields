suppressMessages(here::i_am('R_analysis/analysis/ssa_evi_out_of_sample.R'))
CODE = here::here()
source(glue::glue("{CODE}/R_analysis/init.R"))

plot_ssa_rms = function(oos) {
    
    bar = plot_ssa_rms_bar(oos)
    box = plot_ssa_rms_box(oos)
    
    leg = ggpubr::get_legend(bar)
    return (
        cowplot::plot_grid(box, leg, ncol=1, rel_heights = c(10, 1))
    )
}


plot_ssa_rms_bar = function(out) {
    
    plot_df = out |>
        mutate(rmse = (1 - rmse/baseline)*100) |>
        group_by(crop_station, source) |>
        summarize_all(mean, na.rm=T) |>
        arrange(source, crop_station) |>
        select(-baseline)

    plot_df$source = factor(
        plot_df$source,
        levels = c('cru', 'era5', 'gmfd'),
        labels = c('CRU', 'ERA5-Land', 'GMFD')
    )

    plot_df$crop_station = factor(
        plot_df$crop_station,
        levels = c("all", "only_stations", "only_not_stations"),
        labels = c(
            'All SSA Grid Cells',
            'SSA Grid Cells\nWith Stations',
            'SSA Grid Cells\nWithout Stations'
        )
    )
    
    plot = ggplot(plot_df, aes(x=source, y = rmse, fill=source, alpha=crop_station)) + 
        geom_col(position="dodge", color='black', width=0.7) +
        scale_fill_manual(values = c('#BBBBBB', '#4477AA', '#228833'), guide="none") +
        scale_alpha_manual(values = c(.9,0.5,0.3), name="") +
        scale_y_continuous(expand=c(0,0), breaks=seq(0,21,5), limits=c(0,21)) +
        theme_classic() +
        labs(y='Percent reduction in RMS') +
        theme(
            legend.position='bottom',
            legend.box.background = element_rect(colour = "black"),
            axis.text=element_text(size=18),
            axis.title=element_text(size=18),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line( linewidth=.1, color="grey"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=18),
            strip.background = element_blank(),
            strip.text = element_text(size=18),
            legend.text=element_text(size=18)
            # axis.text.x = element_text(angle=45, vjust=1, hjust=1)
        )     
    return(plot)
}

plot_ssa_rms_box = function(out) {
    
    plot_df = out |>
        mutate(rmse = (1 - rmse/baseline)*100)

    plot_df$source = factor(
        plot_df$source,
        levels = c('cru', 'era5', 'gmfd'),
        labels = c('CRU', 'ERA5-Land', 'GMFD')
    )

    plot_df$crop_station = factor(
        plot_df$crop_station,
        levels = c("all", "only_stations", "only_not_stations"),
        labels = c(
            'All SSA Grid Cells',
            'SSA Grid Cells\nWith Stations',
            'SSA Grid Cells\nWithout Stations'
        )
    )
    means = plot_df |> 
        group_by(source, crop_station) |> 
        summarize(mean_rmse = mean(rmse, na.rm=TRUE), .groups = 'drop')
    plot_df = plot_df |> left_join(means, by=c('source', 'crop_station'))

    plot = ggplot(plot_df, aes(x=source, y = rmse, fill=source, alpha=crop_station)) + 
        geom_boxplot(outlier.shape=NA, coef=0, na.rm=T) +
        geom_point(aes(y=mean_rmse), position=dodge, size=3) +
        scale_fill_manual(values = c('#BBBBBB', '#4477AA', '#228833'), guide="none") +
        scale_alpha_manual(values = c(.9,0.5,0.3), guide="none") +
        # scale_y_continuous(expand=c(0,0), breaks=seq(0,21,5), limits=c(0,21)) +
        theme_classic() +
        labs(y='Percent reduction in RMS') +
        coord_cartesian(ylim = c(0,25)) +
        theme(
            legend.position='bottom',
            legend.box.background = element_rect(colour = "black"),
            axis.text=element_text(size=18),
            axis.title=element_text(size=18),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line( linewidth=.1, color="grey"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=18),
            strip.background = element_blank(),
            strip.text = element_text(size=18),
            legend.text=element_text(size=18)
            # axis.text.x = element_text(angle=45, vjust=1, hjust=1)
        )      
    return(plot)
}