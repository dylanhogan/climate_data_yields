growing_season_min = 3
growing_season_max = 9

year_min = 1950
year_max = 2020

GDD_lower = 10
KDD_lower_list = list(
    corn = list(
        era5=27,
        gmfd=28,
        prism=29,
        era5_sine=27,
        prism_resample=29
    ),
    soy = list(
        era5=28,
        gmfd=29,
        prism=30
    )
)

# if (!is.null(crop) & !is.null(source))
#     KDD_lower = KDD_lower[[crop]][[source]]

temp_lower = -3
temp_upper = 36

g_P = 'prec + prec2'
controls = 'fips + state[t] + state[t2]'
cluster = 'fips + state'