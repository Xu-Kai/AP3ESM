#!/bin/tcsh -f
#ln -sf /disk4/home/lhl/CESM/cesm105_licom/cesm.run/inputdata_licom/domain_licom.permute_cpl7_20120819.nc domain_licom.permute.nc 

setenv link_pathathread /home/export/online1/mdt00/shisuan/swgbcm/licom/data/GB25/file_20250209 
setenv link_path0 /home/export/online1/mdt00/shisuan/swgbcm/licom/data/10km-80/file_20240921 

rm -f Basin.data
rm -f domain_licom.nc
rm -f kmt.da
rm -f temp_initial
rm -f salt_initial
rm -f grid.da
rm -f domain_licom.permute.nc
rm -f tidal_energy.nc
rm -f sss_phc3_monthly.nc
rm -f MODEL_CHLFRC
rm -f dncoef.h1
rm -f depth.da

ln -s $link_pathathread/dncoef.h1_80 dncoef.h1
ln -s $link_pathathread/depth_1_80.da depth.da
ln -s $link_pathathread/temp-20230718-init-map-nn.nc  temp_initial
ln -s $link_pathathread/salt-20230718-init-map-nn.nc  salt_initial
ln -s $link_path0/tidal_energy_licom_tripole_3600X2302X80_omip_db_20240921.nc  tidal_energy.nc
ln -s $link_path0/sss_month_phc3_tripole_3600X2302X80_omip_db_20240921.nc  sss_phc3_monthly.nc
ln -s $link_path0/chl_clim_monthly_tripole_3600X2302X80_omip_db_20240921.nc  MODEL_CHLFRC
ln -s $link_pathathread/LICOM_Hor_Grid_3600X2302X80_tripole_omip_20250209.DATA grid.da
ln -s $link_pathathread/LICOM_Topo_3600X2302X80_tripole_omip_20250209.DATA kmt.da
ln -s $link_pathathread/LICOM_Basin_3600X2302X80_tripole_omip_20250209.DATA Basin.data
ln -s $link_pathathread/domain_licom_tripole_omip_20250209.nc domain_licom.nc

#ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/dataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_10km_80/ocn_001560 .
ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/icedataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_10km_20230718/ocn_001560 .



