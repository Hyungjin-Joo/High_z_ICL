read target

#python 00_pers_mask.py
python 00_astrodrizzle.py

rm *sci1*
rm *single*

ls -1 *flt.fits > flts.lis
python 01_plane_fitting.py

ls -1 *plf.fits > plfs.lis
python 02_WCS_sampling.py

for f in *plf.fits
do
    echo "$f"
    IFS='.'
    read -ra starr<<<"$f"
    sex ${starr[0]}.${starr[1]} -c default.sex -CATALOG_NAME ${starr[0]}.cat -CHECKIMAGE_NAME ${starr[0]}_check.fits
    echo "SExtractor Done"
    echo "$f"
    echo "$f" | python 03_expanding_mask.py
    echo "04"
    echo "$f" | python 04_sky_estimation.py

    rm *check*
    rm *cat
done

#echo "$target" | python 05_astrodrizzle.py
