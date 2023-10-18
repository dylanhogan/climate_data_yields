suppressMessages(here::i_am('R_analysis/figures/climate_data.R'))
CODE = here::here()
if (is.null(init)) source(glue("{CODE}/R_analysis/init.R"))
local_load('analysis/clean_data.R')

generate_climate_data_comparison_df = function(

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