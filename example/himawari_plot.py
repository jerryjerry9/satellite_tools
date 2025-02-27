from ael_satellite_tools.plotting import Himawari
#import matplotlib.pyplot as plt
#import numpy as np

data_path = '/data/C.jerryjerry9/hima_download/himawari_data/'
lon = [90, 180]
lat = [-10, 50]
hima_plot = Himawari(work_path=[],data_path=data_path,plotting_lat_range=lat,plotting_lon_range=lon)

### Set download time period -- 1
year = ['2023'] # Year: from 2015
mon = ['12']    # Month: 01 02 ... 12
day = ['22']    # Day: 01 02 ... 30 31
hour = ['04']   # Hour: 00 01 ... 23
minn = ['00']    # Minute: 00 10 20 30 40 50
# Set band type
AHI_band = [1,2,3,4] # 1, 2, 3 ... 16
geo = ['sun.azm', 'sun.zth','sat.azm', 'sat.zth']

time_list = hima_plot.generate_time_list(year=year,mon=mon,day=day,hour=hour,minn=minn)

file_list, full_path_file_list = hima_plot.generate_data_list(time_list=time_list, AHI_band=AHI_band,geo=geo)

avaiable_time_list, data_issue_list,data_issue_date = hima_plot.check_data(time_list, file_list, full_path_file_list)

file_list, full_path_file_list = hima_plot.generate_data_list(time_list=avaiable_time_list,AHI_band=AHI_band,geo=geo)

output_data_list = []
output_file_list = []
for file_name in full_path_file_list:
    nc_data_list,file_lon_list,file_lat_list,read_file_list = \
        hima_plot.read_nc_file(file_name,missing2nan=True)

    fit_data_array_list, plotting_lon, plotting_lat = \
        hima_plot.fit_resolution(nc_data_list,read_file_list,2,fit_lonlat_output=True)
    output_data_list.append(fit_data_array_list[0])
    output_file_list.append(read_file_list[0])


band_r = output_data_list[2]
band_g = output_data_list[1]
band_b = output_data_list[0]
band_4 = output_data_list[3]
sun_azm = output_data_list[4]
sun_zth = output_data_list[5]
sat_azm = output_data_list[6]
sat_zth = output_data_list[7]

reduce_adjust_angle=78
reduce_adjust = True
band_rr = hima_plot.local_adjustment(band_r,sun_zth,reduce_adjust_angle=reduce_adjust_angle,
                                     reduce_high_zenith_adjust=reduce_adjust)
band_gg = hima_plot.local_adjustment(band_g,sun_zth,reduce_adjust_angle=reduce_adjust_angle,
                                     reduce_high_zenith_adjust=reduce_adjust)
band_bb = hima_plot.local_adjustment(band_b,sun_zth,reduce_adjust_angle=reduce_adjust_angle,
                                     reduce_high_zenith_adjust=reduce_adjust)
band_04 = hima_plot.local_adjustment(band_4,sun_zth,reduce_adjust_angle=reduce_adjust_angle,
                                     reduce_high_zenith_adjust=reduce_adjust)

rs_channel = ['ch3','ch2','ch1'] 
reduce_corr_angle = 78
strength = 1
reduce_corr = True
cor_band_r = hima_plot.rayleigh_correction(avaiable_time_list, sun_azm,sun_zth, 
                                           sat_azm, sat_zth, band_rr,
                                           rs_channel[0], red_band=band_rr,
                                           reduce_corr_angle=reduce_corr_angle, strength=strength,
                                           reduce_rayleigh_corr = reduce_corr)
cor_band_g = hima_plot.rayleigh_correction(avaiable_time_list, sun_azm,sun_zth, 
                                           sat_azm, sat_zth, band_gg,
                                           rs_channel[1], red_band=band_rr,
                                           reduce_corr_angle=reduce_corr_angle, strength=strength,
                                           reduce_rayleigh_corr = reduce_corr)
cor_band_b = hima_plot.rayleigh_correction(avaiable_time_list, sun_azm,sun_zth, 
                                           sat_azm, sat_zth, band_bb,
                                           rs_channel[2], red_band=band_rr,
                                           reduce_corr_angle=reduce_corr_angle, strength=strength,
                                           reduce_rayleigh_corr = reduce_corr)

band_hyg = hima_plot.hybrid_band(cor_band_g ,band_04, data_2_ratio=0.07)

b_rr = cor_band_r
b_gg = band_hyg
b_bb = cor_band_b

min_threshold = [0,0,0]
max_threshold = [100,100,100]
b_r = hima_plot.rescale_value(b_rr,min_threshold[0],max_threshold[0])
b_g = hima_plot.rescale_value(b_gg,min_threshold[1],max_threshold[1])
b_b = hima_plot.rescale_value(b_bb,min_threshold[2],max_threshold[2])

gamma = [1,1,1]
profile_ID = 0
profile_test = None
self_prof= [True,True,True]
enh_r = hima_plot.rgb_enhancement(b_r,gamma=gamma[0],profile_ID=profile_ID, 
                                  self_defined_profile=profile_test, 
                                  self_defined_enhance=self_prof[0])
enh_g = hima_plot.rgb_enhancement(b_g,gamma=gamma[1],profile_ID=profile_ID, 
                                  self_defined_profile=profile_test, 
                                  self_defined_enhance=self_prof[1])
enh_b = hima_plot.rgb_enhancement(b_b,gamma=gamma[2],profile_ID=profile_ID, 
                                  self_defined_profile=profile_test, 
                                  self_defined_enhance=self_prof[2])

rgb_array = hima_plot.rgb_merged(band_red=enh_r,band_green=enh_g,band_blue=enh_b)

rgb_product_name = ['true color']

hima_plot.generate_rgb_figure(rgb_array, plotting_lon, plotting_lat,
                              figure_name=rgb_product_name, time_list=avaiable_time_list,
                              coast_line_color='gold', lonlat_step=4, font_size=22, 
                              prefix='rgb_figure', dpi=300, save_fig=True)


