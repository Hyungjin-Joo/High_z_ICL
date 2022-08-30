from drizzlepac import astrodrizzle

target = input('Target name? ')

astrodrizzle.AstroDrizzle('*plf.fits', skyfile = 'sky_level.txt', output = '%s_F105W'%target, driz_sep_rot = 0, driz_sep_scale = 0.05, driz_sep_outnx = 3558, driz_sep_outny = 3331, combine_type = 'median', final_rot = 0, final_scale = 0.05, final_outnx = 3558, final_outny = 3331, final_units = 'counts')

