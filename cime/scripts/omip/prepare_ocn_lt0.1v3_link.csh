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

setenv DATAPATH_n /home/export/online1/mdt00/shisuan/swgbcm/licom/data/10km-80/file_20250427

ln -sf ${DATAPATH_n}/ss_initial_20230718_for_licom_10km.nc  salt_initial
ln -sf ${DATAPATH_n}/tt_initial_20230718_for_licom_10km.nc  temp_initial

ln -sf ${DATAPATH_n}/uu_initial_20230718_for_licom_10km.nc  uu_initial
ln -sf ${DATAPATH_n}/vv_initial_20230718_for_licom_10km.nc  vv_initial
ln -sf ${DATAPATH_n}/ssh_initial_20230718_for_licom_10km.nc ssh_initial

ln -sf ${DATAPATH_n}/chl_clim_monthly_tripole_3600X2302X80_omip_db_20240921.nc MODEL_CHLFRC
ln -sf ${DATAPATH_n}/sss_month_phc3_tripole_3600X2302X80_omip_db_20240921.nc sss_phc3_monthly.nc
ln -sf ${DATAPATH_n}/tidal_energy_licom_tripole_3600X2302X80_omip_db_20240921.nc tidal_energy.nc

ln -sf ${DATAPATH_n}/depth_L80_20241227.da depth.da
ln -sf ${DATAPATH_n}/LICOM_Hor_Grid_3600X2302X80_tripole_omip_20250427.DATA grid.da
ln -sf ${DATAPATH_n}/LICOM_Topo_3600X2302X80_tripole_omip_20250427.DATA  kmt.da
ln -sf ${DATAPATH_n}/dncoef_L80_20241227.h1 dncoef.h1
ln -sf ${DATAPATH_n}/LICOM_Basin_3600X2302X80_tripole_omip_20250427.DATA Basin.data

ln -sf ${DATAPATH_n}/ice_tripole_omip_20250427.grid icegrid.da
ln -sf ${DATAPATH_n}/ice_tripole_omip_20250427.kmt  icekmt.da

ln -sf ${DATAPATH_n}/domain_licom_tripole_omip_20250427.nc domain_licom_tripole_omip.nc

#ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/icedataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_10km_20230718/ocn_001560 .
ln -sf /home/export/online1/mdt00/shisuan/swgbcm/mxyu/icedataprocess/ocn_init/ocn_init_bilinear/strong/mapping_data_10km_20230718_lt0.1v3/ocn_001560 .



