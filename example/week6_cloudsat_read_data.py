import numpy as np

from ael_satellite_tools.preprocess import CloudSat
from ael_satellite_tools.plotting import Himawari
lat = [-10, 60]
lon = [90, 150]

#cloudsat = CloudSat(work_path=[],lat_range=lat,lon_range=lon)

data_path = '/data/C.jerryjerry9/hima_download/himawari_data'
hima_plot = Himawari(work_path=[],data_path=data_path,plotting_lat_range=lat,plotting_lon_range=lon)

cloudsat = CloudSat(work_path=[],lat_range=lat,lon_range=lon)

print(cloudsat.product_list())

geoprof_file = ['/data/dadm1/obs/CloudSat/hdf-GEOPROF/2016/2016359034029_56694_CS_2B-GEOPROF_GRANULE_P1_R05_E06_F01.hdf']
geoprof_lid_file = ['/data/dadm1/obs/CloudSat/hdf-GEOPROF-LIDAR/2016/2016359034029_56694_CS_2B-GEOPROF-LIDAR_GRANULE_P2_R05_E06_F01.hdf']


Vdata_name,Vdata_len,SD_name = cloudsat.hdf_information(geoprof_file)

for i in range(0,5):
    print(Vdata_name[i],Vdata_len[i])
print(SD_name)

factor, missing_value = cloudsat.cloudsat_var_attribute('Radar_Reflectivity')
print(factor, missing_value)

### read ori data
#ori_lon_list = cloudsat.read_ori_vdata(extracted_geo,'Longitude')
#ori_ref_list = cloudsat.read_ori_sddata(extracted_geo,'Radar_Reflectivity')
### read ori data

rgb_file = ['/data/C.jerryjerry9/hima_download/himawari_data/composite_data/201612240430_true_color_2km_rgb.nc']
rgb_array,plotting_lon_list,plotting_lat_list,output_file_list  = hima_plot.read_rgb_nc_file(rgb_file)

cloudsat.plot_track_w_rgb(geoprof_file[0],rgb_array[0],
                          plotting_lon_list,plotting_lat_list,
                          coast_line_color='gold',
                          lonlat_c='w',lonlat_width=1,lonlat_order=2,
                          lonlat_step=8,font_size=24,
                          figure_title='CloudSat Track',
                          rgb_product='True Color',
                          rgb_time='04:30UTC',
                          prefix='cloudsat_track_w_truecolor',
                          dpi=300, figure_path='cloudsat_fig',save_fig=True)

target_lon_range = [90, 150]
target_lat_range = [15, 30]

#product_name = 'hdf-GEOPROF', 'hdf-GEOPROF-LIDAR'
extracted_geo,extracted_mask_geo = \
cloudsat.sub_domain_check(geoprof_file,target_lon_range,
                          target_lat_range,mask_output=True)
extracted_lid,extracted_mask_lid = \
cloudsat.sub_domain_check(geoprof_lid_file,target_lon_range,
                          target_lat_range,mask_output=True)

extracted_lon = cloudsat.read_vdata(extracted_geo,'Longitude',
                                    extracted_mask_geo,fit_era5_lon=True)
extracted_lat = cloudsat.read_vdata(extracted_geo,'Latitude',
                                    extracted_mask_geo)

extracted_ref = cloudsat.read_sddata(extracted_geo,'Radar_Reflectivity',
                                     extracted_mask_geo)
extracted_CPR_cldmask = cloudsat.read_sddata(extracted_geo,'CPR_Cloud_mask',
                                             extracted_mask_geo)
extracted_height = cloudsat.read_sddata(extracted_geo,'Height',
                                        extracted_mask_geo)

extracted_cfraction = cloudsat.read_sddata(extracted_lid,'CloudFraction',
                                           extracted_mask_lid)

extracted_uncf = cloudsat.read_sddata(extracted_lid,'UncertaintyCF',
                                      extracted_mask_lid)

extracted_cfraction[0][extracted_uncf[0]<=0] = np.nan
extracted_ref[0][extracted_ref[0]==-88.88] = np.nan
extracted_ref[0][extracted_CPR_cldmask[0]<19.5] = np.nan

extracted_hdf_date,extracted_hdf_id,extracted_hdf_granule, \
extracted_region_start_time,extracted_region_end_time \
= cloudsat.read_geometric_info(extracted_geo, extracted_mask_geo)

cloudsat.plot_profile(extracted_ref,extracted_cfraction,'contourf',
                      extracted_lat,5,extracted_height,[0,20],5,
                      5,
                      extracted_hdf_date,extracted_hdf_granule,
                      extracted_region_start_time,extracted_region_end_time,
                      figure_title='CloudSat Profile',
                      prefix='cldsat_prof',figure_path='cloudsat_fig',
                      dpi=300,save_fig=True)

cloudsat.plot_track_w_rgb(geoprof_file[0],rgb_array[0],
                          plotting_lon_list,plotting_lat_list,
                          extracted_lon,extracted_lat, 
                          extracted_hdf_date,extracted_hdf_id,extracted_hdf_granule, 
                          extracted_region_start_time,extracted_region_end_time,
                          extracted_time=False,coast_line_color='gold',
                          lonlat_c='w',lonlat_width=1,lonlat_order=2,
                          lonlat_step=8,font_size=24,
                          figure_title='CloudSat Track',
                          rgb_product='True Color',rgb_time='04:30UTC',
                          prefix='cloudsat_track_w_truecolor',
                          dpi=300, figure_path='cloudsat_fig',save_fig=True)


