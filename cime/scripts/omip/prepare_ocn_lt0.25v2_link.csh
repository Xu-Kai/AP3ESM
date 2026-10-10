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

setenv DATAPATH_n /home/export/online1/mdt00/shisuan/swgbcm/licom/data/25km/file_20250617  
setenv DATAPATH_s /home/export/online1/mdt00/shisuan/swgbcm/licom/data/25km/file_20250617 
setenv grid /home/export/online1/mdt00/shisuan/swgbcm/licom/data/25km/file_20250617 

echo ${DATAPATH_n}


ln -sf ${DATAPATH_s}/salt_jan_phc3_tripole_1440X867X55_omip_db_20250617.nc  salt_initial
ln -sf ${DATAPATH_s}/temp_jan_phc3_tripole_1440X867X55_omip_db_20250617.nc  temp_initial
#ln -sf ${DATAPATH_s}/uu_initial_20000102_10km.nc  uu_initial
#ln -sf ${DATAPATH_s}/vv_initial_20000102_10km.nc  vv_initial
#ln -sf ${DATAPATH_s}/ssh_initial_20000102_10km.nc ssh_initial


ln -sf ${DATAPATH_n}/chl_clim_monthly_tripole_1440X867X55_omip_db_20250617.nc  MODEL_CHLFRC
ln -sf ${DATAPATH_n}/tidal_energy_licom_tripole_1440X867X55_omip_db_20250617.nc tidal_energy.nc

ln -sf ${DATAPATH_n}/sss_month_phc3_tripole_1440X867X55_omip_db_20250617.nc sss_mercator_monthly.nc
ln -sf ${DATAPATH_n}/sst_month_phc3_tripole_1440X867X55_omip_db_20250617.nc sst_ostia_monthly.nc

ln -sf ${grid}/depth.da depth.da
ln -sf ${grid}/dncoef.h1 dncoef.h1
ln -sf ${grid}/LICOM_Hor_Grid_1440X867X55_tripole_omip_20250617.DATA grid.da
ln -sf ${grid}/LICOM_Topo_1440X867X55_tripole_omip_20250617.DATA  kmt.da
ln -sf ${grid}/LICOM_Basin_1440X867X55_tripole_omip_20250617.DATA Basin.data

ln -sf ${grid}/ice_tripole_omip_20250617.grid icegrid.da
ln -sf ${grid}/ice_tripole_omip_20250617.kmt  icekmt.da

ln -sf ${grid}/domain_licom_tripole_omip_20250617.nc domain_licom_tripole_omip.nc

