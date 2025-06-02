import numpy as np
from ael_satellite_tools.preprocess import RO
from ael_satellite_tools.preprocess import CloudSat
from ael_satellite_tools.plotting import Himawari


lat = [-10, 60]
lon = [90, 150]

hima_data_path = '/data/cloud2025/temporary_data'
hima_plot = Himawari(work_path=[],data_path=hima_data_path,plotting_lat_range=lat,plotting_lon_range=lon)
cloudsat = CloudSat(work_path=[],lat_range=lat,lon_range=lon)

#data_path = '/data/dadm1/obs/RO_profile' 
ro = RO(work_path=[],lat_range=lat,lon_range=lon)

ro.ro_information(detail=False)

### Set searching time period 
time_period = ['2019110420','2019110516']

full_path_file_list = ro.generate_list(time_period,
                                       satellite_overlap='both')

extracted_lon_range = [110, 142]
extracted_lat_range = [10, 40]


ro_file_list,ro_lon_list,ro_lat_list \
= ro.sub_domain_check(full_path_file_list,extracted_lon_range,extracted_lat_range,lonlat_list=True)

rgb_file = ['/data/C.jerryjerry9/hima_download/himawari_data/composite_data/201911050440_true_color_2km_rgb.nc']
rgb_array,plotting_lon_list,plotting_lat_list,output_file_list  = hima_plot.read_rgb_nc_file(rgb_file)

geoprof_file = '/data/dadm1/obs/CloudSat/hdf-GEOPROF/2019/2019309034747_72029_CS_2B-GEOPROF_GRANULE_P1_R05_E09_F00.hdf'

cloudsat_lon_range = [90, 150]
cloudsat_lat_range = [20, 30]

extracted_geo,extracted_mask_geo = cloudsat.sub_domain_check(geoprof_file,cloudsat_lon_range,cloudsat_lat_range,mask_output=True)
#extracted_lid,extracted_mask_lid = cloudsat.sub_domain_check(lidar_file,cloudsat_lon_range,cloudsat_lat_range,mask_output=True)
extracted_lon = cloudsat.read_vdata(extracted_geo,'Longitude',extracted_mask_geo,fit_era5_lon=True)
extracted_lat = cloudsat.read_vdata(extracted_geo,'Latitude',extracted_mask_geo)
extracted_hdf_date,extracted_hdf_id,extracted_hdf_granule, \
extracted_region_start_time,extracted_region_end_time \
= cloudsat.read_geometric_info(extracted_geo, extracted_mask_geo)

ro.plot_ro_distribution(ro_file_list,                                              # RO file list
                        rgb_array,plotting_lon_list,plotting_lat_list,             # RGB data
                        geoprof_file,                                              # CloudSat data
                        extracted_lon,extracted_lat,                               # extracted CloudSat data
                        extracted_hdf_date,extracted_hdf_id,extracted_hdf_granule, # extracted CloudSat data
                        extracted_region_start_time,extracted_region_end_time,     # extracted CloudSat data
                        extracted_lon_range=[],extracted_lat_range=[],             # self-defined plotting domain
                        coast_line_color='gold',
                        lonlat_c='w',lonlat_width=0.1,lonlat_order=1,lonlat_step=8,font_size=24,
                        figure_title='RO distribution',rgb_product='True color',rgb_time='04:30UTC',
                        sounding_station=False,prof_num=True,utc_color=True)


temp_prof = ro.read_profile(ro_file_list,'Temp')
moist_prof = ro.read_profile(ro_file_list,'rh')

#vertical_prof = ro.read_profile(ro_file_list,'Pres')
vertical_prof = ro.read_profile(ro_file_list,'MSL_alt')

prof_num = [11,0,7,8]
for i in prof_num:
    print(ro_file_list[i])

ro.plot_ro_profile(ro_file_list,prof_num=4,
                   moist_type='rh',height_type='MSL_alt',
                   height_range=[0,20], height_step=6,
                   prefix='ro_profile', dpi=300, figure_path='ro_fig',save_fig=True)


#sub_plot_num = 2
#prof_num=4
#fig, ax1 = plt.subplots(1, sub_plot_num, figsize=(10, 6.8), tight_layout=True)
#ro.plot_profile_unit(ax1,temp_prof[prof_num],vertical_prof[prof_num],'Temp','RO',
#                     [0,20],6,
#                     sub_plot_num=sub_plot_num,plot_num=1,
#                     linewidth=4,alpha=1)
### comparing data
#ro.plot_profile_unit(ax1,temp_data,vertical_data,'Temp','ERA5',
#                     sub_plot_num=sub_plot_num,plot_num=1,
#                     linewidth=4,alpha=0.8,profile_label=True,axis_setting=False)
