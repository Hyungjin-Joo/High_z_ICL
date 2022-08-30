from drizzlepac import astrodrizzle

astrodrizzle.AstroDrizzle('*flt.fits', output = 'TEST', driz_sep_rot = 0, driz_sep_scale = 0.05, combine_type = 'median',  final_rot = 0, final_scale = 0.05, final_units = 'counts')

