"""
satellite_plotting.py

This module provides functions for reading pre-processed data, 
doing statistic, and plotting

"""

from ..satellite_general import Preparation

class Himawari(Preparation):
    def __init__(self,work_path=None,
                 plotting_lat_range=[-10,50],plotting_lon_range=[90,180],
                 time_period=['20150707'],
                 data_path='/data/cloud2025/temporary_data',data_end_path='day'):
        super().__init__(work_path,plotting_lat_range,plotting_lon_range,time_period)
        from pathlib import Path
        self.plotting_lat_range = plotting_lat_range
        self.plotting_lon_range = plotting_lon_range
        if self.work_path == None or self.work_path == []:
            self.work_path = Path().cwd()
        self.data_path = data_path
        self.data_end_path = data_end_path

    def rgb_composite(self,time_list,rgb_product,ta_resolution=2,plotting_info=True,reduce_local_adjust_angle=78,reduce_high_zenith_adjust=True,reduce_rs_corr_angle=78,reduce_rs_corr_strength=1,reduce_rayleigh_corr = True,hybrid_data2_ratio=0.07,self_gamma=None,profile_ID=0,self_defined_profile=None,self_defined_enhance=None):
        """
        intergate functions for generate rgb-compostie product
        #self.attribute
        #self.generate_data_list
        #self.check_data
        #self.generate_data_list
        #self.read_nc_file 

        produce three band data & geo file

        self.fitting_resolution
        ### optional
        self.local_adjustment
        self.rayleigh_correction 
        self.hybrid_band 
        ###
        #self.rescale_value [threshold_max] [threshold_min] [flag]
        #self.rgb_enhancement [gamma] [flag]
        #self.rgb_merged  
        #return(rgb_array)

        """
        AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,min_threshold,max_threshold,reverse_flag,gamma,self_prof_flag = self.rgb_attribute(rgb_product)       
        file_list, full_path_file_list = self.generate_data_list(time_list=time_list,AHI_band=AHI_band,geo=geo)        
        avaiable_time_list, data_issue_list,data_issue_date = self.check_data(time_list, file_list, full_path_file_list)
#### loop process
        total_rgb_array = []
        plot_lon = []
        plot_lat = []
        for single_time in avaiable_time_list:
            file_list, full_path_file_list = self.generate_data_list(time_list=single_time,AHI_band=AHI_band,geo=geo)
            band_r = self.band_data_generation(band_method[0],full_path_file_list,ta_resolution,latlon_info=False)
            band_g = self.band_data_generation(band_method[1],full_path_file_list,ta_resolution,latlon_info=False)
            band_b,band_lon,band_lat = self.band_data_generation(band_method[2],full_path_file_list,ta_resolution,latlon_info=True)

            for geo_name in geo:
                geo_var = geo_name.replace('.','_')
                var = f'{geo_var}'
                globals()[var] = self.geo_data_generation(geo_var,full_path_file_list,ta_resolution,latlon_info=False)
            if len(band_method) > 3:
                band_04 = self.band_data_generation(band_method[3],full_path_file_list,ta_resolution,latlon_info=False)
                band_04 = self.local_adjustment(band_04,sun_zth,reduce_adjust_angle=reduce_local_adjust_angle,reduce_high_zenith_adjust=reduce_high_zenith_adjust) 
                re_red_band = band_r.copy()
                re_red_band = self.local_adjustment(re_red_band, sun_zth,reduce_adjust_angle=reduce_local_adjust_angle, reduce_high_zenith_adjust=reduce_high_zenith_adjust)

## setting overwrite
            if self_defined_enhance != None:
                self_prof_flag = self_defined_enhance
            if self_gamma != None:
                gamma = self_gamma
        

            bands = [band_r, band_g, band_b]
            bands_functions = [r_functions , g_functions, b_functions]
        
            enh_band = []
            for i in range(0,3):
                band_data = bands[i]
                function_flag = bands_functions[i] 
                if function_flag[0]:
                    band_data = self.local_adjustment(band_data, sun_zth,reduce_adjust_angle=reduce_local_adjust_angle, reduce_high_zenith_adjust=reduce_high_zenith_adjust)
                if function_flag[1]:
                    band_data = self.rayleigh_correction(avaiable_time_list,sun_azm,sun_zth,sat_azm,sat_zth,band_data,rs_channel[i],red_band=re_red_band,reduce_corr_angle=reduce_rs_corr_angle,strength=reduce_rs_corr_strength,reduce_rayleigh_corr = reduce_rayleigh_corr)
                if function_flag[2]:
                    band_data = self.hybrid_band(band_data ,band_04, data_2_ratio=hybrid_data2_ratio)
                if function_flag[3]:
                    band_data = self.rescale_value(band_data,min_threshold[i],max_threshold[i],reverse=reverse_flag[i])
                if function_flag[4]:
                    band_data = self.rgb_enhancement(band_data,gamma=gamma[i],profile_ID=profile_ID,self_defined_profile=self_defined_profile,self_defined_enhance=self_prof_flag[i])
                enh_band.append(band_data)
            
            rgb_array = self.rgb_merged(band_red=enh_band[0],band_green=enh_band[1],band_blue=enh_band[2]) 
            total_rgb_array.append(rgb_array[0])
            plot_lon.append(band_lon)
            plot_lat.append(band_lat)
        if plotting_info:
            return(total_rgb_array,plot_lon,plot_lat,avaiable_time_list)
        else:
            return(total_rgb_array)

    def band_data_generation(self,read_info,file_list,ta_resolution,latlon_info=True):
        file_list = self.check_list(file_list)
        read_method = read_info[-1]
        #print(read_info)
###    method 1: normal read one band; potential bug '_band_xx' exist in the data_path
        if read_method == 0:
            target_file = self.search_target_band_file(read_info[0],file_list)
            nc_data_list,plotting_lon_list,plotting_lat_list,output_file_list = self.read_nc_file(target_file,missing2nan=True)
            fit_data_array_list,local_lon,local_lat = self.fit_resolution(nc_data_list,output_file_list,ta_resolution,fit_lonlat_output=True)
###  
###    method 2: two bands difference
        elif read_method == 1:
            target_file = self.search_target_band_file(read_info[0],file_list)
            nc_data_list_1,plotting_lon_list,plotting_lat_list,output_file_list_1 = self.read_nc_file(target_file,missing2nan=True)

            target_file = self.search_target_band_file(read_info[1],file_list)
            nc_data_list_2,plotting_lon_list,plotting_lat_list,output_file_list = self.read_nc_file(target_file,missing2nan=True)

            nc_data_diff = nc_data_list_1[0] - nc_data_list_2[0]
            fit_data_array_list,local_lon,local_lat = self.fit_resolution(nc_data_diff,output_file_list_1,ta_resolution,fit_lonlat_output=True)
###                             
###    method 3: for band07 refl calculation
        elif read_method == 2:
            target_file = self.search_target_band_file(read_info[0],file_list)
            nc_data_list_07,plotting_lon_list,plotting_lat_list,output_file_list_07 = self.read_nc_file(target_file,missing2nan=True)
            
            target_file_07 = self.search_target_band_file(read_info[1],file_list)
            nc_data_list_13,plotting_lon_list,plotting_lat_list,output_file_list = self.read_nc_file(target_file_07,missing2nan=True) 

            target_file = list(filter(lambda s: 'sun_zth'  in s, file_list))
            nc_data_list_sun_zth,plotting_lon_list,plotting_lat_list,output_file_list = self.read_nc_file(target_file,missing2nan=True)

            nc_data_band07Relf = self.band7Refl(target_file_07,nc_data_list_07,nc_data_list_13,nc_data_list_sun_zth)
            fit_data_array_list,local_lon,local_lat = self.fit_resolution(nc_data_band07Relf,output_file_list_07,ta_resolution,fit_lonlat_output=True)
###
###    method 4: two bands difference but with different resolution
        elif read_method == 3:
            target_file = self.search_target_band_file(read_info[0],file_list)
            nc_data_list_1,plotting_lon_list,plotting_lat_list,output_file_list = self.read_nc_file(target_file,missing2nan=True)
            fit_data_array_list_1 = self.fit_resolution(nc_data_list_1,output_file_list,ta_resolution,fit_lonlat_output=False)
            target_file = self.search_target_band_file(read_info[1],file_list)
            nc_data_list_2,plotting_lon_list,plotting_lat_list,output_file_list = self.read_nc_file(target_file,missing2nan=True)
            fit_data_array_list_2,local_lon,local_lat = self.fit_resolution(nc_data_list_2,output_file_list,ta_resolution,fit_lonlat_output=True)
            fit_data_array_list = [fit_data_array_list_1[0] - fit_data_array_list_2[0]]
###
        if latlon_info:
            return(fit_data_array_list,local_lon,local_lat)
        else:
            return(fit_data_array_list)
   
    def geo_data_generation(self,read_info,file_list,ta_resolution,latlon_info=True):
        file_list = self.check_list(file_list)
        read_info = self.check_list(read_info)
        target_file = list(filter(lambda s: read_info[0]  in s, file_list))
        nc_data_list,plotting_lon_list,plotting_lat_list,output_file_list = self.read_nc_file(target_file,missing2nan=True)
        fit_data_array_list,local_lon,local_lat = self.fit_resolution(nc_data_list,output_file_list,ta_resolution,fit_lonlat_output=True)
        if latlon_info:
            return(fit_data_array_list,local_lon,local_lat)
        else:
            return(fit_data_array_list)

    def generate_data_list(self,time_list,product_name=[],AHI_band=[],geo=[],band4km=[],separate_ori_4km=False):
        ## self-defined or given product name for auto-generating
        time_list = self.check_list(time_list)
        band_table, AHI_band_num_table, band_resolution, geo_name_table, nc_data_type_table_4km = self.string_info(nc_plotting_list=True)
        if len(product_name) > 0:
            AHI_band,geo = self.rgb_attribute(product_name,band_info_only=True)
            band4km=[]
        file_list = []
        full_path_file_list = []
        for file_date in time_list:
            #path_year = file_date[0:4]
            #path_mon = file_date[4:6]
            #nc_data_path = self.data_path + '/sub_domain_data/'+path_year+'/'+path_mon+''
            date_path = self.date_path_generate(file_date)
            nc_data_path = self.data_path + '/sub_domain_data/'+date_path+''
            file_name = []
            if len(AHI_band) > 0:
                for b_num in AHI_band:
                    if b_num > 16:
                        break
                    band_num = str(b_num)
                    if b_num < 10:
                         band_num = '0' + band_num
                    file_name = file_date+'_band_' + band_num +'.nc'
                    file_list.append(file_name)
                    full_path_file_list.append(nc_data_path+'/'+file_name)
                    pos_idx = AHI_band_num_table.index(band_num)
                    band_type_num = band_table[pos_idx]
                    band_type = band_type_num[0:3]
                    if len(band4km) > 0:
                        for band_var in band4km:
                            if (band_var == 'rad' or band_var == 'tbb') and band_type == 'tir':
                                file_name = file_date+'_4km_band_' + band_num +'_'+band_var+'.nc'
                                file_list.append(file_name)
                                full_path_file_list.append(nc_data_path+'/'+file_name)
                            if (band_var == 'rad' or band_var == 'rfc' or band_var == 'rfy') and band_type != 'tir':
                                file_name = file_date+'_4km_band_' + band_num +'_'+band_var+'.nc'
                                file_list.append(file_name)
                                full_path_file_list.append(nc_data_path+'/'+file_name)
            if len(geo) > 0:
                for geo_name in geo:
                    pos_idx = geo_name_table.index(geo_name)
                    file_name = file_date + '_4km_'+nc_data_type_table_4km[pos_idx]+'.nc'
                    file_list.append(file_name)
                    full_path_file_list.append(nc_data_path+'/'+file_name) 
        return(file_list, full_path_file_list)
             
    def check_data(self, time_list, file_list, full_path_file_list):
        import glob
        time_list = self.check_list(time_list)
        file_list = self.check_list(file_list)
        full_path_file_list = self.check_list(full_path_file_list)
        data_issue_list = []
        data_issue_date = []
        for i in range(0,len(full_path_file_list)):
        ## nc exist or not
           checking_file = sorted(glob.glob(full_path_file_list[i]))
           if len(checking_file) < 1: 
               data_issue_list.append(file_list[i])
               file_date_list = self.nc_name_info(file_list[i],date_info=True)
               if file_date_list[0] in time_list:
                   data_issue_date.append(file_date_list[0])
                   time_list = list(filter(lambda x: file_date_list[0] not in x, time_list))

        ## nc exist, check lat lon range
           else:
              band_resol_list = self.nc_name_info(file_list[i],resol_info=True)
              domain_flag,nc_lon_range,nc_lat_range = self.read_nc_boundary(checking_file,self.plotting_lon_range,self.plotting_lat_range,band_resol_list[0])
              if domain_flag == False:
                  data_issue_list.append(file_list[i])
                  file_date_list = self.nc_name_info(file_list[i],date_info=True)
                  if file_date_list[0] in time_list:
                      data_issue_date.append(file_date_list[0])
                      time_list = list(filter(lambda x: file_date_list[0] not in x, time_list))
        ## generate available time list for following plotting process
        ## generate data issue list for pre-process to check data status
        return(time_list, data_issue_list,data_issue_date)

    def read_nc_file(self,full_path_file_list,missing2nan=True):
        import netCDF4 as nc
        import numpy as np
        full_path_file_list = self.check_list(full_path_file_list)
        nc_data_list = []
        plotting_lon_list = []
        plotting_lat_list = []
        output_file_list = []
        for file_name in full_path_file_list:
            print('read file:'+file_name+'')
            band_resol_list = self.nc_name_info(file_name,resol_info=True)
        ## read nc var and latlon info
            local_lon,local_lat = self.lonlat_index(self.plotting_lon_range, self.plotting_lat_range,band_resol_list[0],lonlat_only=True)
            nc_re = nc.Dataset(file_name, 'r',  format='NETCDF4_CLASSIC')
            var_names = list(nc_re.variables.keys())
            nclon = list(nc_re.variables['lon'])
            nclat = list(nc_re.variables['lat'])
            lon_idx_s = nclon.index(local_lon[0])
            lon_idx_e = nclon.index(local_lon[-1])
            lat_idx_s = nclat.index(local_lat[0])
            lat_idx_e = nclat.index(local_lat[-1])
            plotting_lon = np.array(nclon[lon_idx_s:lon_idx_e+1])
            plotting_lat = np.array(nclat[lat_idx_s:lat_idx_e+1])
            nc_var_name = var_names[-1]
            nc_data = np.array(nc_re.variables[nc_var_name][0,lat_idx_s:lat_idx_e+1,lon_idx_s:lon_idx_e+1])
            if missing2nan:
                 missing_value = getattr(nc_re.variables[nc_var_name], "missing_value", None)
                 nc_data[nc_data==missing_value] = np.nan
            nc_re.close()
            nc_data_list.append(nc_data)
            plotting_lon_list.append(plotting_lon)
            plotting_lat_list.append(plotting_lat)
            output_file_list.append(file_name)
        return(nc_data_list,plotting_lon_list,plotting_lat_list,output_file_list) 

    def read_rgb_nc_file(self,file_list,missing2nan=True):
        import netCDF4 as nc
        import numpy as np
        file_list = self.check_list(file_list)
        nc_data_list = []
        plotting_lon_list = []
        plotting_lat_list = []
        output_file_list = []
        for file_name in file_list:

            print('read file:'+file_name+'')
            split_name = file_name.split('/')[-1]
            split_parts = split_name.replace('.','_').split('_')
            #rgb_var = split_parts[-2].split('')
            resol_name = list(filter(lambda s: 'km'  in s, split_parts))
            len_pos = 0 - len(resol_name[0])
              
            resol = resol_name[0][len_pos:-2]
            zero_num = 0
            for test_0 in resol:
                if int(test_0) == 0:
                   zero_num = zero_num +1
            resol = int(resol)/10**zero_num
            resolution = int(12000/resol)
            band_resol_list = int(resolution)
        ## read nc var and latlon info
            local_lon,local_lat = self.lonlat_index(self.plotting_lon_range, self.plotting_lat_range,band_resol_list,lonlat_only=True)
            nc_re = nc.Dataset(file_name, 'r',  format='NETCDF4_CLASSIC')
            var_names = list(nc_re.variables.keys())
            nclon = list(nc_re.variables['lon'])
            nclat = list(nc_re.variables['lat'])
            lon_idx_s = nclon.index(local_lon[0])
            lon_idx_e = nclon.index(local_lon[-1])
            lat_idx_s = nclat.index(local_lat[0])
            lat_idx_e = nclat.index(local_lat[-1])
            plotting_lon = np.array(nclon[lon_idx_s:lon_idx_e+1])
            plotting_lat = np.array(nclat[lat_idx_s:lat_idx_e+1])

            data_array = []
            for i in range(3,len(var_names)):
                nc_data = np.array(nc_re.variables[var_names[i]][0,lat_idx_s:lat_idx_e+1,lon_idx_s:lon_idx_e+1])
                if missing2nan:
                     missing_value = getattr(nc_re.variables[var_names[i]], "missing_value", None)
                     nc_data[nc_data==missing_value] = np.nan
                data_array.append(nc_data)
            nc_re.close()
#####       
            if len(data_array) == 3: 
                merged_array = self.rgb_merged(band_red=data_array[0],band_green=data_array[1],band_blue=data_array[2])
                nc_data_list.append(merged_array[0])
            else:
                nc_data_list.append(data_array)
##### 
            plotting_lon_list.append(plotting_lon)
            plotting_lat_list.append(plotting_lat)
            output_file_list.append(file_name)
        return(nc_data_list,plotting_lon_list,plotting_lat_list,output_file_list)

    """
    image adjustment function
    self.fit_resolution
    self.local_adjustment
    self.rayleigh_correction
    self.hybrid_band
    self.rescale_value
    self.rgb_enhancement
    """
    def fit_resolution(self,nc_data_list,file_list,ta_resolution,fit_lonlat_output=False):
        import cv2
        import numpy as np
        file_list = self.check_list(file_list)
        nc_data_list = self.check_data_list(nc_data_list)
#        plotting_lon_list = self.check_list(plotting_lon_list)
#        plotting_lat_list = self.check_list(plotting_lat_list)
        fit_data_array_list = []
        print('fit data into '+str(ta_resolution)+'km resolution')
        for i in range(0,len(nc_data_list)):
            resol = []
#            if len(plotting_lon_list) > 0:
#                plotting_lon = plotting_lon_list[i]
#                if len(plotting_lon) > 1:
#                    resol = round((plotting_lon[1] - plotting_lon[0])*1000)/10
#            if len(plotting_lat_list) > 0 and len(resol) < 1:
#                plotting_lat = plotting_lat_list[i]
#                if len(plotting_lat) > 1:
#                    resol = round((plotting_lat[1] - plotting_lat[0])*1000)/10
            if len(file_list) > 0:
                band_resol_list = self.nc_name_info(file_list[i],resol_info=True)
                resol = 12000/band_resol_list[0]
            scale_factor = resol/ta_resolution
###    fitting resolution
            data_array = nc_data_list[i]
            array_size = data_array.shape
            ## some 4km data should not be applied to resize method
            if scale_factor != 1:
                print('fit resolution')
                ## angle data should apply different resize method  
                split_name = file_list[i].split('/')[-1]
                split_parts = split_name.replace('_', '.').split('.')
                if '4km' in split_parts:
                    if 'sun' in split_parts or 'sat' in split_parts:
                        angles = np.deg2rad(data_array)  
                        sin_vals = np.sin(angles)
                        cos_vals = np.cos(angles)
                        resized_sin = cv2.resize(sin_vals, (int(array_size[1]*scale_factor), int(array_size[0]*scale_factor)), interpolation=cv2.INTER_LINEAR)
                        resized_cos = cv2.resize(cos_vals, (int(array_size[1]*scale_factor), int(array_size[0]*scale_factor)), interpolation=cv2.INTER_LINEAR)
                        fit_data_array = np.rad2deg(np.arctan2(resized_sin, resized_cos)) % 360
                    elif 'rfy' in split_parts or  'tbb' in split_parts:  
                        fit_data_array = cv2.resize(data_array, (int(array_size[1]*scale_factor), int(array_size[0]*scale_factor)), interpolation=cv2.INTER_LINEAR)
                    else:
                        print('currently not provide the resize data')   
                        print('output ori-resolution data')   
                        fit_data_array = data_array       
                else:
                    fit_data_array = cv2.resize(data_array, (int(array_size[1]*scale_factor), int(array_size[0]*scale_factor)), interpolation=cv2.INTER_LINEAR)
                fit_data_array_list.append(fit_data_array)
            else:
                print('output ori-resolution data') 
                fit_data_array_list.append(data_array)
###    output data
        if fit_lonlat_output:
            local_lon,local_lat = self.lonlat_index(self.plotting_lon_range, self.plotting_lat_range,resolution=ta_resolution,lonlat_only=True) 
            return(fit_data_array_list,local_lon,local_lat)
        else:
            return(fit_data_array_list)

    def local_adjustment(self,data_array_list,sun_zenith_angle_list,reduce_adjust_angle=78,reduce_high_zenith_adjust=True):
        import numpy as np
        data_array_list = self.check_list(data_array_list)
        sun_zenith_angle_list = self.check_list(sun_zenith_angle_list) 
        adjusted_data = []
        print('Adjust local brightness')
        for i in range(0,len(data_array_list)):
            #if np.max(data_array_list[i])>121 and force_adjust=False:
            #    print('Input data may not be reflectivity...')
            #    print('Output non-adjusted data')
            #    adjusted_data.append(data_array_list[i])
            #else:
                adjust_angle = sun_zenith_angle_list[i].copy()
                if reduce_high_zenith_adjust:
                    adjust_angle[adjust_angle > reduce_adjust_angle] = reduce_adjust_angle
                sun_zenith_rad = np.radians(adjust_angle)
                local_adjust = np.cos(sun_zenith_rad)
                adjusted_data.append(data_array_list[i]/local_adjust)
        return(adjusted_data)
        
    def rayleigh_correction(self,time_list,sun_azm,sun_zth,sat_azm,sat_zth,band_list,ch,red_band=None,reduce_corr_angle=78,strength=1,reduce_rayleigh_corr = True):
        from pyspectral.rayleigh import Rayleigh
        import numpy as np
        if red_band == None:
            N = len(band_list)
            red_band = [None] * N
        band_list = self.check_list(band_list)
        time_list = self.check_list(time_list)
        red_band = self.check_list(red_band)
        sun_azm = self.check_list(sun_azm)
        sun_zth = self.check_list(sun_zth)
        sat_azm = self.check_list(sat_azm)
        sat_zth = self.check_list(sat_zth)
        corrected_band_data = []

        for i in range(0,len(band_list)):
            band_date = int(time_list[i])
            if band_date > 202212130449:
                sat_name = 'Himawari-9'
            elif band_date > 201802130249 and band_date < 201802140710:
                sat_name = 'Himawari-9'
            else:
                sat_name = 'Himawari-8'
        ### for Rayleigh correction: remove the bluish haze
            
            az_diff = sun_azm[i] - sat_azm[i]
            hima = Rayleigh(sat_name, 'ahi')
#####################
            if ch == 'ch1':
               print('Correct blue band data')
            elif ch == 'ch2':
               print('Correct green band data')
            elif ch == 'ch3':
               print('Correct red band data')
            refl_cor_band = hima.get_reflectance(sun_zth[i], sat_zth[i], az_diff, ch, red_band[i])
            if reduce_rayleigh_corr:
                thresh_zen = reduce_corr_angle
                maxzen = np.nanmax(sun_zth[i])
                strength = 1
                refl_cor_band_reduce = hima.reduce_rayleigh_highzenith(sun_zth[i], refl_cor_band, thresh_zen, maxzen,strength)
                corrected_band_data.append(band_list[i] - refl_cor_band_reduce)
            else:
                corrected_band_data.append(band_list[i] - refl_cor_band)
#######################
        return(corrected_band_data)

    def hybrid_band(self,band_data_1,band_data_2,data_2_ratio=0.07):
        data_1_ratio = 1-0.07
        band_data_1 = self.check_list(band_data_1)
        band_data_2 = self.check_list(band_data_2)
        hybrid_data = []
        print('Combine two band data with second bata ratio = '+str(data_2_ratio)+'')
        for i in range(0,len(band_data_1)):
            hybrid_data.append(band_data_1[i]*data_1_ratio + band_data_2[i]*data_2_ratio)
        return(hybrid_data)
 
    def rescale_value(self,data_array_list,min_threshold,max_threshold,reverse=False):
        data_array_list = self.check_list(data_array_list)
        min_threshold = self.check_list(min_threshold)
        max_threshold = self.check_list(max_threshold)
        rescaled_data = []
        rescale_min = min_threshold[0]
        rescale_max = max_threshold[0]
        rescale_range = rescale_max - rescale_min
        rescaled_data = []
        print('Ignore the value outside the range of '+str(rescale_min)+' ~ '+str(rescale_max)+'')
        print('Data are rescaled to the range of 0 ~ 1 depending on scale range')
        for i in range(0,len(data_array_list)):

            single_data = data_array_list[i]
            single_data[single_data < rescale_min] = rescale_min
            single_data[single_data > rescale_max] = rescale_max
            if reverse == False:
                single_data = (single_data-rescale_min)/rescale_range
            else:
                single_data = abs(single_data-rescale_max)/rescale_range
            rescaled_data.append(single_data)
        return(rescaled_data)
            
    def rgb_enhancement(self,data_array_list,gamma=None,profile_ID=None,self_defined_profile=None,self_defined_enhance=False):
        import numpy as np
        data_array_list = self.check_list(data_array_list)
        enh_data = []
#        nc_data[nc_data==missing_value] = np.nan

        if self_defined_enhance:
            print('Apply self-defined enhancement method')
            enh_profile = self.self_defined_enh_profile(profile_ID)
            if self_defined_profile != None:
                print('User-defined enhancement profile')
                enh_profile = self_defined_profile
            for i in range(0,len(data_array_list)):
                for_nan = data_array_list[i].copy()
                single_data = self.enh_truecolor_profile(data_array_list[i],enh_profile)
                single_data = np.where(np.isnan(for_nan), np.nan, single_data)
                enh_data.append(single_data)
        else:
            print('Gamma enhancement value is set to '+str(gamma)+'')
            for i in range(0,len(data_array_list)):
                for_nan = data_array_list[i].copy()
                single_data = (data_array_list[i])**(1/gamma)
                single_data = np.where(np.isnan(for_nan), np.nan, single_data)
                enh_data.append(single_data)
        return(enh_data)

    """
    plotting or output data
    self.rgb_merged
    self.generate_rgb_figure
    self.generate_rgb_nc_file
    """
    def rgb_merged(self,band_red=None,band_green=None,band_blue=None):
        import numpy as np
        band_red = self.check_list(band_red)
        band_green = self.check_list(band_green)
        band_blue = self.check_list(band_blue)

        ref_band = next((k for k in (band_red, band_green, band_blue) if k is not None), None)
        total_rgb_array = []
        for i in range(0,len(ref_band)):
            array_size = ref_band[i].shape
            rgb_array=np.empty((array_size[0],array_size[1],3))
            rgb_array[:,:,0] = band_red[i]
            rgb_array[:,:,1] = band_green[i]
            rgb_array[:,:,2] = band_blue[i]
            total_rgb_array.append(rgb_array)
        return(total_rgb_array) 

   
    def generate_rgb_figure(self,rgb_array, plotting_lon, plotting_lat,figure_name=[],time_list=[],coast_line_color='gold',lonlat_step=4,font_size=24,prefix='rgb_figure',dpi=300, figure_path=[], save_fig=True):
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
        from pathlib import Path
        rgb_array = self.check_list(rgb_array)
        figure_name = self.check_list(figure_name)
        time_list = self.check_list(time_list)
        plotting_lon = self.check_list(plotting_lon)
        plotting_lat = self.check_list(plotting_lon)
########        
        lon_s = self.plotting_lon_range[0]
        lon_e = self.plotting_lon_range[1]
        lat_s = self.plotting_lat_range[0]
        lat_e = self.plotting_lat_range[1]
        lon_range = lon_e - lon_s
        lat_range = lat_e - lat_s
        lon_step = np.floor((lon_range)/(lonlat_step-1))
        lat_step = np.floor((lat_range)/(lonlat_step-1))
######## 
        for i in range(0,len(rgb_array)):
            plot_array = np.where(np.isnan(rgb_array[i]), 0, rgb_array[i])

            lon2,lat2=np.meshgrid(plotting_lon[i], plotting_lat[i])
            ## regular imshow plotting
            fig,ax1 = plt.subplots(1,1,figsize=(13,13.5*(lat_range/lon_range)),tight_layout=True)
            m = Basemap(llcrnrlon=lon_s, urcrnrlon=lon_e, llcrnrlat=lat_s, urcrnrlat=lat_e,resolution='i')
            m.drawparallels(np.arange(lat_s, lat_e+1, lat_step), labels=[1, 0, 0, 0], linewidth=0.01, color='k', fontsize=font_size-2,zorder=1)
            m.drawmeridians(np.arange(lon_s, lon_e+1, lon_step), labels=[0, 0, 0, 1], linewidth=0.01, color='k', fontsize=font_size-2,zorder=1)
            m.drawcoastlines(linewidth=1.,color=coast_line_color,zorder=3)
            
            m.imshow(plot_array,zorder=2)
            title_date = ['']
            if len(time_list)>0:
                date = time_list[i]
                print(date)
                title_date = [''+date[0:4]+'/'+date[4:6]+'/'+date[6:8]+' '+date[8:10]+':'+date[10:12]+'UTC']
            rgb_name = ['']
            output_name = ['']
            if len(figure_name)>0:
                rgb_name = [figure_name[0].replace('_',' ').title()]
                output_name = ['_'+ figure_name[0].replace(' ','_').lower() +'_']
            plt.title(''+rgb_name[0]+' '+title_date[0]+'',fontsize=font_size)
            if save_fig:
               if figure_path == None or figure_path == []:
                   figure_path = ''+str(self.work_path)+'/'
               else:
                   Path(figure_path).mkdir(parents=True, exist_ok=True)
                   figure_path = ''+str(figure_path)+'/'
               plt.savefig(''+figure_path+''+prefix+''+output_name[0]+''+date+'.png',dpi=dpi)

    def generate_rgb_nc_file(self,time_list,plotting_lon,plotting_lat,ta_resolution,product_name,domain_name=[],rgb_array=[],r_band=[],g_band=[],b_band=[],nc_output=True):
        import os
        import pickle
        import numpy as np
        from netCDF4 import Dataset
        from numpy import dtype
        time_list = self.check_list(time_list)
        resolu_name = str(ta_resolution).replace('.','')
        plotting_lon = self.check_list(plotting_lon)
        plotting_lat = self.check_list(plotting_lat)
        #file_name = self.check_list(file_name)
        rgb_array = self.check_list(rgb_array)
        r_band = self.check_list(r_band)
        g_band = self.check_list(g_band)
        b_band = self.check_list(b_band)
        product = product_name[0].replace(' ','_')
        product = product.lower()
        print('Generating nc file')
        for i in range(0,len(time_list)):
            nc_data_path = self.data_path + '/composite_data'
            os.makedirs(nc_data_path, exist_ok=True)
            data_date = time_list[i]
            band_var = ''
            r_var = True
            g_var = True
            b_var = True
###    band var naming rule
            if len(rgb_array)>0:
                out_r, array_size = self.convert_nan2value(rgb_array[i][:,:,0],-999,shape_return=True)
                out_g = self.convert_nan2value(rgb_array[i][:,:,1],-999)
                out_b = self.convert_nan2value(rgb_array[i][:,:,2],-999)
                band_var = band_var + 'rgb'
            else:
                if len(r_band)>0:
                    out_r, array_size = self.convert_nan2value(r_band[i],-999,shape_return=True)
                    band_var = band_var + 'r'
                else:
                    r_var = False
                if len(g_band)>0:
                    out_g, array_size = self.convert_nan2value(g_band[i],-999,shape_return=True)
                    band_var = band_var + 'g'
                else:
                    g_var = False
                if len(b_band)>0:
                    out_b, array_size = self.convert_nan2value(b_band[i],-999,shape_return=True)
                    band_var = band_var + 'b'
                else:
                    b_var = False
            if domain_name == []:
                output_name = ''+data_date+'_'+product+'_'+str(resolu_name)+'km_'+band_var+''
            else:
                domain = domain_name[0].replace(' ','_')
                output_name = ''+data_date+'_'+product+'_'+str(resolu_name)+'km_'+band_var+'_'+str(domain)+''
###
            local_lon = plotting_lon[i]
            local_lat = plotting_lat[i]

            if nc_output:
                tin = np.double(np.arange(0,2,1))
                ncout = Dataset(''+nc_data_path+'/'+output_name+'.nc', 'w', format='NETCDF4')
                ncout.createDimension('time', None)  # unlimited
                ncout.createDimension('lat', array_size[0])
                ncout.createDimension('lon', array_size[1])
               #create time axis
                time = ncout.createVariable('time', dtype('double').char, ('time',))
                time.long_name = 'UTC time'
                time.units = 'days since '+data_date[0:4]+'-'+data_date[4:6]+'-'+data_date[6:8]+' '+data_date[8:10]+':'+data_date[10:12]+':00'
                time.calendar = 'standard'
                time.axis = 'T'
               # create latitude axis
                lat = ncout.createVariable('lat', dtype('float').char, ('lat'))
                lat.standard_name = 'latitude'
                lat.long_name = 'latitude'
                lat.units = 'degrees_north'
                lat.axis = 'Y'
               # create longitude axis
                lon = ncout.createVariable('lon', dtype('float').char, ('lon'))
                lon.standard_name = 'longitude'
                lon.long_name = 'longitude'
                lon.units = 'degrees_east'
                lon.axis = 'X'

               # create variable array
                if r_var:
                    dout1 = ncout.createVariable('red_band', dtype('float').char, ('time', 'lat', 'lon'))
                    dout1.long_name = 'Rescaled RGB data: red band'
                    dout1.units = '0~1'
                    dout1.missing_value = -999

                if g_var:
                    dout2 = ncout.createVariable('green_band', dtype('float').char, ('time', 'lat', 'lon'))
                    dout2.long_name = 'Rescaled RGB data: green band'
                    dout2.units = '0~1'
                    dout2.missing_value = -999
                if b_var:
                    dout3 = ncout.createVariable('blue_band', dtype('float').char, ('time', 'lat', 'lon'))
                    dout3.long_name = 'Rescaled RGB data: blue band'
                    dout3.units = '0~1'
                    dout3.missing_value = -999

               # copy axis from original dataset
                time[:] = tin[0]
                lon[:] = local_lon[:]
                lat[:] = local_lat[:]
                if r_var:
                    output1 = np.empty((1,array_size[0],array_size[1]))
                    output1[0,:,:] = out_r
                if g_var:
                    output2 = np.empty((1,array_size[0],array_size[1]))
                    output2[0,:,:] = out_g
                if b_var:
                    output3 = np.empty((1,array_size[0],array_size[1]))
                    output3[0,:,:] = out_b
                if r_var:
                    dout1[:] = output1[:]
                if g_var:
                    dout2[:] = output2[:]
                if b_var:
                    dout3[:] = output3[:]
                ncout.close()
                print(''+output_name+'.nc generated')
            else:
               # output .pkl file
                with open(''+output_name+'.pkl', 'wb') as f:
                    pickle.dump(output_array, f)
                print(''+output_name+'.pkl generated')

    """
    preparation sub function
    """
    def band7Refl(self,band07_file_list,band07_list,band13_list,sun_zenith_angle_list):
        import numpy as np
        c39,a,b,c,c1,c2,v = self.band7Refl_constant()
        band07_file_list = self.check_list(band07_file_list)
        band07_list = self.check_list(band07_list)
        band13_list = self.check_list(band13_list)
        sun_zenith_angle_list = self.check_list(sun_zenith_angle_list)
        band07_relf = []
        for i in range(0,len(band07_list)):
            date = self.nc_name_info(band07_file_list[i],date_info=True)
            year = int(date[0][0:4])
            mon = int(date[0][4:6])
            day = int(date[0][6:8])
            jd = int(self.calculate_julian_date(year,mon,day) - 0.5)
            band07 = band07_list[i]
            band13 = band13_list[i]
            sun_zenith_angle = sun_zenith_angle_list[i]

            sun_zth_2km = self.fit_resolution(sun_zenith_angle,'4km_sun_zth',2)
            ESD = 1-0.0167*(np.cos((2*np.pi*(jd-3))/365))
            TOARAD = (c39/ESD**2)*np.cos(2*np.pi*(sun_zth_2km[0]/360))
            Rtot =  (c1*(v**3))/(np.exp((c2*v)/(a + b*band07 + c*band07**2))-1)
            Rtherm = (c1*(v**3))/(np.exp((c2*v)/(a + b*band13 + c*band13**2))-1)
            relf = 100*((Rtot-Rtherm)/(TOARAD-Rtherm))
            band07_relf.append(relf)
        return(band07_relf)

    def search_target_band_file(self,band_num,file_list):
        file_list = self.check_list(file_list)
        read_num = str(band_num)
        if band_num <10:
            read_num = '0' + read_num
        read_band = '_band_' + read_num
        target_file = list(filter(lambda s: read_band  in s, file_list))
        return(target_file)

    def nc_name_info(self,file_list,date_info=False,resol_info=False):
        file_list = self.check_list(file_list)
        band_table, AHI_band_num_table, band_resolution, geo_name_table, nc_data_type_table_4km = self.string_info(nc_plotting_list=True)
        band_resol_list = []
        file_date_list = []
        for file_name in file_list:
            split_name = file_name.split('/')[-1]
            split_parts = split_name.split('_')
            if date_info:
                file_date_list.append(split_parts[0])
            if resol_info:
                if '4km' in split_parts:
                    band_resol_list.append(3000)
                else:
                    band_num_idx = split_parts.index('band') + 1
                    band_num = split_parts[band_num_idx][0:2]
                    pos_idx = AHI_band_num_table.index(band_num)
                    band_resol_list.append(band_resolution[pos_idx])
        if date_info:
            return(file_date_list)
        if resol_info:
            return(band_resol_list)

    def convert_nan2value(self,data_array,value,shape_return=False):
        import numpy as np
        output = np.nan_to_num(data_array, nan=value)
        array_size = output.shape
        if shape_return:
            return(output,array_size)
        else:
            return(output)

    def enh_truecolor_profile(self,band_data,enh_profile):
        import numpy as np
        band_data = band_data*255
 
        data_array_shape = band_data.shape
        lat_size = data_array_shape[0]
        lon_size = data_array_shape[1]
        enh_band = np.zeros((lat_size,lon_size))
        enh_band[:,:] = 255
        x = enh_profile[0]
        y = enh_profile[1]
        for i in range(0,4):
            x1 = x[i]
            x2 = x[i+1]
            y1 = y[i]
            y2 = y[i+1]
            m = (y2-y1)/(x2-x1)
            b = y2-(m*x2)
            mask1 = np.zeros((lat_size,lon_size))
            mask1[band_data>=x1] = 1
            mask2 = np.zeros((lat_size,lon_size))
            mask2[band_data<x2] = 1
            new_mask = mask1*mask2
            enh_band[new_mask>0] = band_data[new_mask>0]*m + b
        enh_band = enh_band/255
        return(enh_band)

    def enh_truecolor_profile_3D(self,band_data,data_array_shape):
        import numpy as np
        band_data = band_data*255
        lat_size = data_array_shape[0]
        lon_size = data_array_shape[1]
        enh_band = np.zeros((lat_size,lon_size,3))
        enh_band[:,:,:] = 255
        profile_ID = 0
        enh_profile = self.self_defined_enh_profile(profile_ID)
        x = enh_profile[0]
        y = enh_profile[1]
        for i in range(0,4):
            x1 = x[i]
            x2 = x[i+1]
            y1 = y[i]
            y2 = y[i+1]
            m = (y2-y1)/(x2-x1)
            b = y2-(m*x2)
            mask1 = np.zeros((lat_size,lon_size,3))
            mask1[band_data>=x1] = 1
            mask2 = np.zeros((lat_size,lon_size,3))
            mask2[band_data<x2] = 1
            new_mask = mask1*mask2
            enh_band[new_mask>0] = band_data[new_mask>0]*m + b
        enh_band = enh_band/255
        return(enh_band)
    
    def self_defined_enh_profile(self,profile_ID):
###    self-defined profile 
        profile_modis = [[0, 25,  55, 100, 255],[0, 90, 140, 175, 255]]
###    intergate all profiles for ID targeting
        total_profile = [profile_modis]
        return(total_profile[profile_ID])

    def band7Refl_constant(self):
        c39 =  4.8077
        a = 0.4793907798197780
        b = 0.999234381214647
        c = 1.85684785537253*(10**(-7))
        c1 = 1.19104*(10**(-5))
        c2 = 1.43878
        v = 2575.767
        return(c39,a,b,c,c1,c2,v)

    """
    RGB composite product
    self.rgb_attribute # calling product info
    rgb attribute information
    """
    def rgb_attribute(self,rgb_product_name,band_info_only=False):
        rgb_product_name = self.check_list(rgb_product_name)
        for rgb_product in rgb_product_name:
            rgb_product = rgb_product.replace(' ', '_')
            rgb_product = rgb_product.lower()
            print(rgb_product)
            method = getattr(self, rgb_product)
            product_attr = method()
            AHI_band = product_attr[0]
            geo = product_attr[1]
        if band_info_only:
           return(AHI_band, geo)
        else:
           return(product_attr)
    
    def true_color(self):
        AHI_band = [1,2,3,4]
        band_method = [[3,0],[2,0],[1,0],[4,0]]
        geo = ['sun.azm', 'sun.zth','sat.azm', 'sat.zth']
        r_functions = [True,True,False,True,True]
        g_functions = [True,True,True,True,True]
        b_functions = [True,True,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [0,0,0]
        thres_max = [100,100,100]
        thres_rev = [False,False,False]
        gamma = [1,1,1]
        self_prof = [True,True,True]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def cloud_phase_distinction(self):
        AHI_band = [3,5,13]
        band_method = [[13,0],[3,0],[5,0]]
        geo = ['sun.zth']
        r_functions = [False,False,False,True,True]
        g_functions = [True,False,False,True,True]
        b_functions = [True,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [219.6,0,1]
        thres_max = [280.7,85,50]
        thres_rev = [True,False,False]
        gamma = [1,1,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def day_deep_clouds(self):
        AHI_band = [3,8,13]
        band_method = [[13,8,1],[3,0],[13,0]]
        geo = ['sun.zth']
        r_functions = [False,False,False,True,True]
        g_functions = [True,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-5,70,243.6]
        thres_max = [35,100,292.6]
        thres_rev = [True,False,False]
        gamma = [1,1,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)
              
    def day_snow_fog(self):
        AHI_band = [4,5,7,13]
        band_method = [[4,0],[5,0],[7,13,2]]
        geo = ['sun.zth']
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [0,0,2]
        thres_max = [102,68,45]
        thres_rev = [False,False,False]
        gamma = [1.6,1.7,1.95]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def microphysics_24hr_b14(self):
        AHI_band = [11,13,14,15]
        band_method = [[13,15,1],[14,11,1],[13,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-3,-0.4,248.6]
        thres_max = [7.5,6.1,303.2]
        thres_rev = [True,False,False]
        gamma = [1,1.1,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def microphysics_24hr_b13(self):
        AHI_band = [11,13,15]
        band_method = [[13,15,1],[13,11,1],[13,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-3,0.8,248.6]
        thres_max = [7.5,5.8,303.2]
        thres_rev = [True,False,False]
        gamma = [1,1.3,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def night_microphysics(self):
        AHI_band = [7,13,15]
        band_method = [[13,15,1],[7,13,1],[13,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-3,-7,243.7]
        thres_max = [7.5,2.9,293.2]
        thres_rev = [True,True,False]
        gamma = [1,1,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def day_microphysics_warm(self):
        AHI_band = [4,7,13]
        band_method = [[4,0],[7,13,2],[13,0]]
        geo = ['sun.zth']
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [0,2,203.5]
        thres_max = [102,82,303.2]
        thres_rev = [False,False,False]
        gamma = [0.95,2.6,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def day_microphysics_cold(self):
        AHI_band = [4,7,13]
        band_method = [[4,0],[7,13,2],[13,0]]
        geo = ['sun.zth']
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [0,2,203.5]
        thres_max = [102,38,303.2]
        thres_rev = [False,False,False]
        gamma = [0.95,1.8,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def day_cloud_phase(self):
        AHI_band = [1,5,6]
        band_method = [[5,0],[6,0],[1,0]]
        geo = ['sun.zth']
        r_functions = [True,False,False,True,True]
        g_functions = [True,False,False,True,True]
        b_functions = [True,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [0,0,0]
        thres_max = [50,50,100]
        thres_rev = [False,False,False]
        gamma = [1,1,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def natural_color(self):
        AHI_band = [3,4,5]
        band_method = [[5,0],[4,0],[3,0]]
        geo = ['sun.zth']
        r_functions = [True,False,False,True,True]
        g_functions = [True,False,False,True,True]
        b_functions = [True,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [0,0,0]
        thres_max = [99,102,100]
        thres_rev = [False,False,False]
        gamma = [1,0.95,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def air_mass(self):
        AHI_band = [8,10,12,13]
        band_method = [[10,8,1],[13,12,1],[8,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [0,-4.3,208]
        thres_max = [25.8,41.5,242.6]
        thres_rev = [True,True,True]
        gamma = [1,1,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def differential_water_vapor(self):
        AHI_band = [8,10]
        band_method = [[10,8,1],[10,0],[8,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-3,213.2,208.5]
        thres_max = [30,278.2,243.9]
        thres_rev = [True,True,True]
        gamma = [3.5,2.5,2.5]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def simple_water_vapor(self):
        AHI_band = [8,10,13]
        band_method = [[13,0],[8,0],[10,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [202.3,214.7,245.1]
        thres_max = [279,242.7,261]
        thres_rev = [True,True,True]
        gamma = [10,5.5,5.5]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def day_convective_storms(self):
        AHI_band = [3,5,7,8,10,13]
        band_method = [[10,8,1],[7,13,1],[5,3,3]]
        geo = ['sun.zth']
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [True,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-5,-1,-80]
        thres_max = [36,61,26]
        thres_rev = [True,False,False]
        gamma = [1,0.5,0.95]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def so2_b13(self):
        AHI_band = [9,10,11,13]
        band_method = [[9,10,1],[13,11,1],[13,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-6,-1.6,243.6]
        thres_max = [5,4.9,303.2]
        thres_rev = [False,False,False]
        gamma = [1,1.2,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def so2_b14(self):
        AHI_band = [9,10,11,13,14]
        band_method = [[9,10,1],[14,11,1],[13,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-6,-5.9,243.6]
        thres_max = [5,5.1,303.2]
        thres_rev = [False,False,False]
        gamma = [1,0.85,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def ash_b13(self):
        AHI_band = [11,13,15]
        band_method = [[13,15,1],[13,11,1],[13,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-3,-1.6,243.6]
        thres_max = [7.5,4.9,303.2]
        thres_rev = [True,False,False]
        gamma = [1,1.2,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def ash_b14(self):
        AHI_band = [11,13,14,15]
        band_method = [[13,15,1],[14,11,1],[13,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-3,-5.9,243.6]
        thres_max = [7.5,5.1,303.2]
        thres_rev = [True,False,False]
        gamma = [1,0.85,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def dust_b13(self):
        AHI_band = [11,13,15]
        band_method = [[13,15,1],[13,11,1],[13,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-3,0.9,261.5]
        thres_max = [7.5,12.5,289.2]
        thres_rev = [True,False,False]
        gamma = [1,2.5,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def dust_b14(self):
        AHI_band = [11,13,14,15]
        band_method = [[13,15,1],[14,11,1],[13,0]]
        geo = []
        r_functions = [False,False,False,True,True]
        g_functions = [False,False,False,True,True]
        b_functions = [False,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [-3,-0.5,261.5]
        thres_max = [7.5,15,289.2]
        thres_rev = [True,False,False]
        gamma = [1,2.2,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def natural_fire_color(self):
        AHI_band = [3,4,6]
        band_method = [[6,0],[4,0],[3,0]]
        geo = ['sun.zth']
        r_functions = [True,False,False,True,True]
        g_functions = [True,False,False,True,True]
        b_functions = [True,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [0,0,0]
        thres_max = [100,100,100]
        thres_rev = [False,False,False]
        gamma = [1,1,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def fire_temperature(self):
        AHI_band = [5,6,7]
        band_method = [[7,0],[6,0],[5,0]]
        geo = ['sun.zth']
        r_functions = [False,False,False,True,True]
        g_functions = [True,False,False,True,True]
        b_functions = [True,False,False,True,True]
        rs_channel = ['ch3','ch2','ch1']
        thres_min = [273,0,0]
        thres_max = [350,50,50]
        thres_rev = [False,False,False]
        gamma = [1,1,1]
        self_prof = [False,False,False]
        return(AHI_band,geo,band_method,r_functions,g_functions,b_functions,rs_channel,thres_min,thres_max,thres_rev,gamma,self_prof)

    def rgb_composite_name(self):
        print('RGB composite product name')
        print('-------------------------------------------------------')
        print('Day time only:')
        print('True color; Natural color; Day deep clouds; Day snow fog')
        print('Cloud phase distinction; Day cloud phase; Day convective storms')
        print('Day microphysics warm; Day microphysics cold')
        print('-------------------------------------------------------')
        print('Night time only:')
        print('Night microphysics')
        print('-------------------------------------------------------')
        print('All time:')
        print('Differential water vapor; Simple water vapor; Air mass')
        print('Microphysics 24hr b13; Microphysics 24hr b14')
        print('-------------------------------------------------------')
        print('Target driven')
        print('Natural fire color; Fire temperature')
        print('SO2 b13; SO2 b14; Ash b13; Ash b14')
        print('Dust b13; Dust b14')
        print('-------------------------------------------------------')
        print(' \'warm\', \'cold\', \'b13\', and \'b14\' represent different formulas of product')
