suppressMessages(here::i_am('R_analysis/figures/climate_impacts.R'))
CODE = here::here()
if (is.null(init)) source(glue::glue("{CODE}/R_analysis/init.R"))
local_load('analysis/climate.R')
local_load('analysis/clean_data.R')

generate_climate_impacts_plot_df = function(
    year_min, year_max,
    growing_season_min, growing_season_max,
    GDD_lower, KDD_lower_list,
    temp_lower, temp_upper,
    aggregate=TRUE, t_list=c(1,2,3,4),
    seed=123, regression_args=list(),
    num_iterations=1000, num_cores=30,
    spec_list = c('piecewise', 'bins', 'poly', 'tavg'),
    sources = c("prism", "era5", "gmfd"),
    crops = c("corn", "soy")
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
                temp_lower, temp_upper
            )
            
            # Run climate impacts calculation
            sp_list = list()
            for (spec in spec_list) {
            
                climpacts = run_climate_projection(
                    df, spec, 
                    GDD_lower, KDD_lower, 
                    temp_lower, temp_upper,
                    aggregate, t_list,
                    num_iterations, num_cores, 
                    seed, regression_args
                ) |>
                mutate(crop=!!crop, source=!!source)
                
                df_list = append(df_list, list(climpacts))
            }
        }
    }
                    
    clim_out = bind_rows(df_list)
    return(clim_out)
}
            
plot_climate_impacts = function(plot_df) {
    
    plot_df = plot_df |> 
        select(-iteration) |> 
        as.data.table() |>
        melt(
            id.vars = c('spec', 'crop', 'source'), 
            variable.name='temp', 
            value.name='impact') |>
        group_by(spec, crop, source, temp) |>
        summarize(
            across(
                where(is.numeric), 
                list(
                    mean = ~ mean(.x), 
                    ub = ~ quantile(.x, 0.95),
                    lb = ~ quantile(.x, 0.05)
                )
            ),
            .groups = 'drop'
        ) |>
        ungroup()

    plot_df$source = factor(
        plot_df$source,
        levels = c('prism', 'era5', 'gmfd'),
        labels = c('PRISM', 'ERA5-Land', 'GMFD')
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
    plot_df$temp = factor(
        plot_df$temp,
        levels = c('impact1', 'impact2', 'impact3', 'impact4'),
        labels = c('1°C', '2°C', '3°C', '4°C')
    )

    
    plot = ggplot(plot_df) +
        geom_linerange(
            aes(x = spec, ymin = impact_lb*100, ymax = impact_ub*100, color = source),
            lwd = 1.2, 
            position=position_dodge(width=0.5)
        ) +
        geom_point(
            aes(x = spec, y=impact_mean*100, color=source), 
            position=position_dodge(width=0.5)
        ) +
        facet_wrap(
            ~source, 
            strip.position = "bottom", 
            scales = "free_x"
        ) + 
        facet_grid2(
            vars(crop), 
            vars(temp)
        ) +
        labs(color='Climate dataset', y='Percent change in yield', x='Functional form') +
        scale_color_manual(values=c('#CC3311', '#0077BB', '#009988')) +
        theme_classic() +
        theme(
            legend.position='bottom',
            legend.box.background = element_rect(colour = "black"),
            text = element_text(size = 14),
            plot.title = element_text(hjust = 0.5),
            panel.grid.major.y = element_line( size=.1, color="gray"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle=45, vjust=1, hjust=1),
            strip.text.x = element_text(size = 12),
            strip.text.y = element_text(size = 12),
        )
    
    return(plot)
    
}

plot_climate_map = function(plot_df, shp_dir) {
    
    shp = st_read(shp_dir)

    shp = shp |> 
        filter(!STATEFP %in% c(
            "02","15", "60",
            "66","72","78", "69")) |>
        mutate(fips=as.numeric(GEOID)) |> 
        select(fips, geometry)

    # This loop is annoyingly slow, but necessary to get
    # polygons w/o crops included for each crop-source combo
    plot_gdf = data.frame()
    for (crop in unique(plot_df$crop)) {
        for (source in unique(plot_df$source)) {
            plot_gdf = bind_rows(
                plot_gdf,
                plot_df |>
                    filter(crop==!!crop, source==!!source) |>
                    right_join(shp, by='fips') |>
                    mutate(
                        centroid = st_centroid(geometry), 
                        longitude = st_coordinates(centroid)[, 1],
                        latitude = st_coordinates(centroid)[, 2],
                        crop = !!crop,
                        source = !!source
                    ) |>
                    select(-centroid, -geometry) |>
                    filter(longitude > -100)
            )
        }
    }
    
    plot_df = shp |> right_join(plot_gdf, by='fips')
    
    plot_df$source = factor(
        plot_df$source,
        levels = c('prism', 'era5', 'gmfd'),
        labels = c('PRISM', 'ERA5-Land', 'GMFD')
    )
    plot_df$crop = factor(
        plot_df$crop,
        levels = c('corn', 'soy'),
        labels = c('Corn', 'Soybeans')
    )

    # Topcode (bottomcode?) results. This only affects a handful of 
    # regions and makes plots more consistent and readable
    plot_df = plot_df |> mutate(impact = ifelse(impact2 < -0.5, -0.5, impact2))
    
    plot = ggplot() +
        geom_sf(data=plot_df, aes(fill = impact*100)) +
        scale_fill_gradient2(
            high="#ca0020", low='#0571b0', mid='#f7f7f7', 
            limits= c(-50, 50), breaks=seq(-50,50,10)
        ) + 
        facet_grid2(vars(source),vars(crop)) +
        labs(fill="Percent change\n      in yield") +
        guides(fill = guide_colorbar(barwidth = unit(12, "cm"))) +
        theme_classic() + 
        theme(
            legend.direction = "horizontal",
            legend.position='bottom',
            axis.title = element_blank(),  # Remove axis titles
            axis.text = element_blank(),   # Remove axis text labels
            axis.ticks = element_blank(),   # Remove axis ticks
            panel.grid = element_blank(),   # Remove grid lines
            axis.line = element_blank(),
            text = element_text(size = 16))
    
    return(plot)
}