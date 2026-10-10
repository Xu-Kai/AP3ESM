#!/bin/tcsh -f

setenv DATAPATH_s /home/export/online1/mdt00/shisuan/swgbcm/licom/data/GB25/licom_2km_250728
setenv restore /home/export/online1/mdt00/shisuan/swgbcm/licom/data/GB25/licom_2km_250728/data
setenv inifield /home/export/online1/mdt00/shisuan/swgbcm/licom/data/GB25/licom_2km_250728/ini_field


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


ln -sf ${inifield}/ss_initial_20230723_for_licom_2km.nc  salt_initial
ln -sf ${inifield}/tt_initial_20230723_for_licom_2km.nc  temp_initial
ln -sf ${inifield}/uu_initial_20230723_for_licom_2km.nc  uu_initial
ln -sf ${inifield}/vv_initial_20230723_for_licom_2km.nc  vv_initial
ln -sf ${inifield}/ss_initial_20230723_for_licom_2km.nc  ssh_initial

ln -sf ${restore}/chl_clim_monthly_tripole_18000X11511X80_omip_db_20240922.nc MODEL_CHLFRC
ln -sf ${restore}/SSS_mercator_monthly_2023_for_licom_2km.nc sss_mercator_monthly.nc
ln -sf ${restore}/tidal_energy_licom_tripole_18000X11511X80_omip_db_20240922.nc tidal_energy.nc

ln -sf ${DATAPATH_s}/depth_L80_20241227.da depth.da
ln -sf ${DATAPATH_s}/LICOM_Hor_Grid_18000X11511X80_tripole_omip_20250505.DATA grid.da
ln -sf ${DATAPATH_s}/LICOM_Topo_18000X11511X80_tripole_omip_20250505.DATA  kmt.da
ln -sf ${DATAPATH_s}/dncoef_L80_20241227.h1 dncoef.h1
ln -sf ${DATAPATH_s}/LICOM_Basin_18000X11511X80_tripole_omip_20250505.DATA Basin.data

ln -sf ${DATAPATH_s}/ice_tripole_omip_20250505.grid icegrid.da
ln -sf ${DATAPATH_s}/ice_tripole_omip_20250505.kmt  icekmt.da

ln -sf ${DATAPATH_s}/domain_licom_tripole_omip_20250505.nc domain_licom_tripole_omip.nc

ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/gbdataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_2km_0723/ocn_028800 .



