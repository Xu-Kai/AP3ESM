#!/bin/tcsh -f

setenv DATAPATH_s /home/export/online1/mdt00/shisuan/swgbcm/licom/data/GB25/file_20250220

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

ln -sf ${DATAPATH_s}/salt-20230718-init-map-nn.nc  salt_initial
ln -sf ${DATAPATH_s}/temp-20230718-init-map-nn.nc  temp_initial
ln -sf ${DATAPATH_s}/chl_clim_monthly_tripole_10800X6907X80_omip_db_20250220.nc MODEL_CHLFRC
ln -sf ${DATAPATH_s}/sss_month_phc3_tripole_10800X6907X80_omip_db_20250220.nc sss_phc3_monthly.nc
ln -sf ${DATAPATH_s}/tidal_energy_licom_tripole_10800X6907X80_omip_db_20250220.nc tidal_energy.nc
ln -sf ${DATAPATH_s}/depth_1_80.da depth.da
ln -sf ${DATAPATH_s}/LICOM_Hor_Grid_10800X6907X80_tripole_omip_20250220.DATA grid.da
ln -sf ${DATAPATH_s}/LICOM_Topo_10800X6907X80_tripole_omip_20250220.DATA  kmt.da
ln -sf ${DATAPATH_s}/dncoef.h1_80 dncoef.h1
ln -sf ${DATAPATH_s}/LICOM_Basin_10800X6907X80_tripole_omip_20250220.DATA Basin.data
ln -sf ${DATAPATH_s}/domain_licom_tripole_omip_20250220.nc domain_licom_tripole_omip.nc

ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/icedataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_3km_80/ocn_010440 .
#ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/dataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_2km_80/ocn_028800 .



