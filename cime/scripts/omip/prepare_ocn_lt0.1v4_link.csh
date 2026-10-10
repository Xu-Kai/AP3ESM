#!/bin/tcsh -f
#ln -sf /disk4/home/lhl/CESM/cesm105_licom/cesm.run/inputdata_licom/domain_licom.permute_cpl7_20120819.nc domain_licom.permute.nc 

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

setenv DATAPATH_n /home/export/online1/mdt00/shisuan/swgbcm/licom/data/10km-80/xjw_20250827/phc3_month 
setenv DATAPATH_s /home/export/online1/mdt00/shisuan/swgbcm/licom/data/10km-80/xjw_20250827/ini_field/20000102 
setenv grid /home/export/online1/mdt00/shisuan/swgbcm/licom/data/10km-80/xjw_20250827/readin_files_20250812


ln -sf ${DATAPATH_s}/ss_20000102_for_swlicom_forecast_10km.nc  salt_initial
ln -sf ${DATAPATH_s}/tt_20000102_for_swlicom_forecast_10km.nc  temp_initial
ln -sf ${DATAPATH_s}/uu_20000102_for_swlicom_forecast_10km.nc  uu_initial
ln -sf ${DATAPATH_s}/vv_20000102_for_swlicom_forecast_10km.nc  vv_initial
ln -sf ${DATAPATH_s}/ssh_20000102_for_swlicom_forecast_10km.nc ssh_initial

ln -sf ${DATAPATH_n}/chl_clim_monthly_10km.nc MODEL_CHLFRC
ln -sf ${DATAPATH_n}/tidal_energy_licom_tripole_3600X2302X80_omip_db_20240921.nc tidal_energy.nc

ln -sf ${DATAPATH_n}/sss_month_phc3_14.nc sss_mercator_monthly.nc
ln -sf ${DATAPATH_n}/sst_month_phc3_14.nc sst_ostia_monthly.nc


ln -sf ${grid}/depth_L80_20241227.da depth.da
ln -sf ${grid}/dncoef_L80_20241227.h1 dncoef.h1
ln -sf ${grid}/LICOM_Hor_Grid_3600X2302X80_tripole_omip_20250809.DATA grid.da
ln -sf ${grid}/LICOM_Topo_3600X2302X80_tripole_omip_20250809.DATA  kmt.da
ln -sf ${grid}/LICOM_Basin_3600X2302X80_tripole_omip_20250809.DATA Basin.data

ln -sf ${grid}/ice_tripole_omip_20250809.grid icegrid.da
ln -sf ${grid}/ice_tripole_omip_20250809.kmt  icekmt.da

ln -sf ${grid}/domain_licom_tripole_omip_20250809.nc domain_licom_tripole_omip.nc

ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/newgitdataprocess/ocn_init/ocn_init_bilinear/strong/data_10km_lt0.1v4_20000102/ocn_001560 .
ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/newgitdataprocess/ocn_init/ocn_init_bilinear/strong/data_10km_lt0.1v4_20000102/ocn_002808 .


