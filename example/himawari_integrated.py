from ael_satellite_tools.preprocess import Himawari as Himawari
from ael_satellite_tools.plotting import Himawari as Hima_plot

lat = [-10, 50]
lon = [90, 180]
data_path = '/data/C.jerryjerry9/hima_download/himawari_data'

himawari = Himawari(work_path=[],data_path=data_path,lat_range=lat,lon_range=lon)
hima_plot = Hima_plot(work_path=[],data_path=data_path,plotting_lat_range=lat,plotting_lon_range=lon)

### Set download time period -- 1
year = ['2018'] # Year: from 2015
mon = ['01']    # Month: 01 02 ... 12
day = ['10']    # Day: 01 02 ... 30 31
hour = ['03']   # Hour: 00 01 ... 23
minn = ['10']    # Minute: 00 10 20 30 40 50
### Set download time period -- 2
#time_period = ['202312220400','202312240400']
#time_delta = ['days=1','hours=0','minutes=0']

### Set band type
### AHI band number
#AHI_band = [1,2,3,4] # 1, 2, 3 ... 16
### CEReS band type & number
#band = ['VIS']  # band: VIS TIR SIR EXT
#band_num = [1]      # band number: 1 2 ... 10
### geo data
#geo = ['sun.azm', 'sun.zth','sat.azm', 'sat.zth']

hima_plot.rgb_composite_name()
rgb_product_name = ['Day microphysics warm']
AHI_band, geo = hima_plot.rgb_attribute(rgb_product_name,band_info_only=True)

file_list = []
zip_file_list = []
ftp_path_file_list = []
for band_name in AHI_band:
    band, band_num = himawari.band_name_convert(band_name)
    tem_file_list, tem_zip_file_list, tem_ftp_path_file_list = \
    himawari.generate_list(time_period=[],time_delta=[],year=year,mon=mon,day=day,
                           hour=hour,minn=minn,band=band,band_num=band_num)
    file_list.extend(tem_file_list)
    zip_file_list.extend(tem_zip_file_list)
    ftp_path_file_list.extend(tem_ftp_path_file_list)
tem_file_list, tem_zip_file_list, tem_ftp_path_file_list = \
himawari.generate_list(time_period=[],time_delta=[],year=year,mon=mon,day=day,
                       hour=hour,minn=minn,geo=geo,band4km=[])
file_list.extend(tem_file_list)
zip_file_list.extend(tem_zip_file_list)
ftp_path_file_list.extend(tem_ftp_path_file_list)


himawari.pre_process(ftp_path_file_list,remove_list_flag=True)

AHI_band,geo,band_method,\
r_functions,g_functions,b_functions,rs_channel,\
min_threshold,max_threshold,reverse_flag,\
gamma,self_enh_flag = hima_plot.rgb_attribute(rgb_product_name)

time_list = hima_plot.generate_time_list(time_period=[],time_delta=[],
                                         year=year,mon=mon,day=day,hour=hour,minn=minn)

reduce_local_adjust_angle=78
reduce_high_zenith_adjust = True
reduce_rs_corr_angle=78
reduce_rs_corr_strength=1
reduce_rayleigh_corr = True
ta_resolution = 2
self_gamma=None
profile_ID=0
self_defined_profile=None
self_defined_enhance=None

rgb_array, plotting_lon, plotting_lat, avaiable_time_list = \
hima_plot.rgb_composite(time_list, rgb_product_name, ta_resolution=2, plotting_info=True,
                        reduce_local_adjust_angle=78, reduce_high_zenith_adjust=True,
                        reduce_rs_corr_angle=78, reduce_rs_corr_strength=1, 
                        reduce_rayleigh_corr = True,
                        hybrid_data2_ratio=0.07,
                        self_gamma=None,profile_ID=0,
                        self_defined_profile=None,self_defined_enhance=None)

hima_plot.generate_rgb_nc_file(avaiable_time_list,plotting_lon,plotting_lat,ta_resolution,
                               rgb_product_name,rgb_array=rgb_array,
                               domain_name=[],nc_output=True)

hima_plot.generate_rgb_figure(rgb_array, plotting_lon, plotting_lat,
                              figure_name=rgb_product_name, time_list=avaiable_time_list,
                              coast_line_color='gold', lonlat_step=4, font_size=22, 
                              prefix='rgb_figure', dpi=100, 
                              figure_path='himawari_figure', save_fig=True)

#rgb,plotting_lon_list,plotting_lat_list,output_file_list =\
#hima_plot.read_rgb_nc_file(full_path_rgb_file_list)
