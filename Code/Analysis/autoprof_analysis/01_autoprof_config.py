target = 'XMM1229'

ap_process_mode = 'image'

ap_image_file = 'detection_samp.fits'
ap_name = '%s_detection'%target
ap_pixscale = 0.05
ap_guess_center = {'x': 912., 'y': 744.}
ap_mask_file = '%s_seg_BCG_banned_2.0_add_samp.fits'%target
ap_set_psf = 5
ap_doplot = True
ap_isoclip = True
ap_truncate_evaluation = True
