00_autoprof_cut.py: crop the drizzled image to run the AutoProf
01_autoprof_config.py: configuration file for AutoProf
02_autoprof_sfb.py: estimate the surface brightness profile
03_autoprof_plot.py: convert the surface brightness in [counts/arcs2] to [mag/arcs2] and plot the profile
04_alt_autoprof_pmn.py: do multi-Sérsic profile fitting with pymultinest. U need to correct the code to modify the number of components.
06_autoprof_pmn_fraction_icl.py: based on pymultinest result, estimate the ICL flux and uncertainties.