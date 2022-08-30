read target
echo "$target" | python 01_detection_file.py
echo "$target" | python 02_rms_file.py
for f in *sci.fits
do
    IFS='.'
    read -ra starr<<<"$f"

    sex detection.fits,${starr[0]}.fits -c default.sex -CATALOG_NAME ${starr[0]}.cat -CHECKIMAGE_NAME ${target}_seg.fits -BACK_SIZE 128 -DETECT_MINAREA 3 -DETECT_THRESH 1.5 -WEIGHT_IMAGE map_weight.fits,${starr[0]}_rms.fits -WEIGHT_TYPE MAP_WEIGHT,MAP_RMS -WEIGHT_THRESH 0.01,1.e+6
done
echo "SExtractor Done"
echo "$target" | python 03_expanding_mask.py
echo "Expanding mask Done"
echo "$target" | python 04_sky_dsky.py
