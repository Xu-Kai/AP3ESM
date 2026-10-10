#!/bin/tcsh -f

setenv DATAPATH_s /home/export/online1/mdt00/shisuan/swgbcm/licom/data/GB25/file_20250222
setenv path2 /home/export/online1/mdt00/shisuan/swgbcm/licom/data/2km-80/file_20250221

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

ln -sf ${DATAPATH_s}/salt_20230718_remap-nn.nc  salt_initial
ln -sf ${DATAPATH_s}/temp-20230718_remap-nn.nc  temp_initial
ln -sf ${DATAPATH_s}/chl_clim_monthly_tripole_36000X23023X80_omip_db_20250222.nc MODEL_CHLFRC
ln -sf ${DATAPATH_s}/sss_month_phc3_tripole_36000X23023X80_omip_db_20250222.nc sss_phc3_monthly.nc
ln -sf ${DATAPATH_s}/tidal_energy_licom_tripole_36000X23023X80_omip_db_20250222.nc tidal_energy.nc
ln -sf ${path2}/depth_1_80.da depth.da
ln -sf ${DATAPATH_s}/LICOM_Hor_Grid_36000X23023X80_tripole_omip_20250222.DATA grid.da
ln -sf ${DATAPATH_s}/LICOM_Topo_36000X23023X80_tripole_omip_20250222.DATA  kmt.da
ln -sf ${path2}/dncoef.h1_80 dncoef.h1
ln -sf ${DATAPATH_s}/LICOM_Basin_36000X23023X80_tripole_omip_20250222.DATA Basin.data
ln -sf ${DATAPATH_s}/domain_licom_tripole_omip_20250222.nc domain_licom_tripole_omip.nc

#ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/dataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_2km_80/ocn_028800 .
ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/icedataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_1km_80/ocn_057600 .
#ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/icedataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_1km_80/ocn_115200 .



