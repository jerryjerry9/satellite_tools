"""
satellite_tools.py

This module provides functions downloading & pre-processing satellite products
pre-processing may include extract sub-domain, re-grid, composite ...

planning satellite product: Himawari, GPM-IMERG, CloudSat, FS7, ASCAT ...

"""

class Preparation:
    def __init__(self,work_path,lat_range,lon_range,time_period):
        self.work_path = work_path
        self.lat_range = lat_range
        self.lon_range = lon_range
        self.time_period = time_period


class Data_group:
    def __init__(self,band_data):
        self.band_data = band_data


class Himawari(Preparation):
    def __init__(self,work_path=None,lat_range=[-10,50],lon_range=[90,180],time_period=['20150707'],data_path='/data/dadm1/obs/Himawari'):
        super().__init__(work_path,lat_range,lon_range,time_period)
        self.data_path = data_path

    def test_part(self):
        data_path = self.data_path
        print(data_path)

    def generate_list(self,time_period=[],time_delta=[],year=[],mon=[],day=[],hour=[],minn=[],band=[],band_num=[],band4km=[],geo=[],separate_ori_4km=False):
        """
        Integrate functions for generating list for download
        self.generate_time_list
        self.generate_data_list
        """
###    generate time list from continious time range or specific time 
        time_list = self.generate_time_list(time_period=time_period,time_delta=time_delta,year=year,mon=mon,day=day,hour=hour,minn=minn)
###    generate data list for downloading data
        file_list, zip_file_list, full_path_file_list = self.generate_data_list(time_list=time_list,band=band,band_num=band_num,geo=geo,band4km=band4km,separate_ori_4km=separate_ori_4km)
        return (file_list, zip_file_list, full_path_file_list)

    def pre_process(self, full_path_file_list, remove_list_flag=True, download_flag=True):
        """
        Integrate all functions for pre-processing
        download data, sub-domain extraction, and generate nc files
        self.check_exist_sub_domain_file
        self.download
        self.unzip
        self.read_binary
        self.convert
        self.sub_domain_extract
        self.generate_nc
        """
        full_path_file_list = self.check_exist_sub_domain_file(full_path_file_list,remove_list_flag=remove_list_flag)

        for file_name in full_path_file_list:

            output_file_list = self.download(file_name, download_flag=download_flag)
            output_file_list, data_array = self.unzip(output_file_list)  
            # generate_binary(output_file_list,data_array)
            output_file_list, np_binary_data = self.read_binary(output_file_list, data_array)

            band_tbb = self.convert(output_file_list, np_binary_data)
            sub_band_tbb, sub_lon, sub_lat = self.sub_domain_extract(band_tbb)
            self.generate_nc(output_file_list, sub_band_tbb, sub_lon, sub_lat)

    def generate_data_list(self,time_list=[],band=[],band_num=[],band4km=[],geo=[],separate_ori_4km=False):
###  prepare vars
        self.__himawari_FTP = 'ftp://hmwr829gr.cr.chiba-u.ac.jp/gridded/FD/V20190123'
        self.band = band
        self.band_num = band_num
        self.band4km = band4km
        self.geo  = geo
###  generating download flag
        if len(band) > 0 and len(band_num) > 0 and separate_ori_4km == False:      
            down_ori = True
        else:
            down_ori = False
        if len(band4km) > 0:
            down_4km = True
        else:
            down_4km = False
        if len(geo) > 0:
            geo_info = True
        else:
            geo_info = False
###  time setting
###  return date list        
        #time_list = self.generate_time_list(time_period=time_period,time_delta=time_delta,year=year,mon=mon,day=day,hour=hour,minn=minn)
###  generating file list for download from FTP           
        file_list = []
        zip_file_list = [] 
        full_path_file_list = []
###    time period
        for single_time in time_list:
            YYYY = single_time[0:4] 
            MM = single_time[4:6]
            DD = single_time[6:8]
            HH = single_time[8:10]
            MN = single_time[10:12]
###    band type
            for CHN in band:
                        for NUM in band_num:
                            if CHN == 'VIS' and NUM > 3:
                                break
                            elif CHN == 'SIR' and NUM > 2:
                                break
                            elif CHN == 'EXT' and NUM > 1:
                                break
                            if NUM < 10:
                                num = '0'+str(NUM)  
###    ori-resolution band data
                            if down_ori:
                                file_name = [''+YYYY+''+MM+''+DD+''+HH+''+MN+'.'+CHN.lower()+'.'+num+'.fld.geoss']
                                download_file = [file_name[0] + '.bz2']
                                full_path_file = [''+self.__himawari_FTP+'/'+YYYY+MM+'/'+CHN+'/'+download_file[0]+'']
                                file_list.append(file_name[0])
                                zip_file_list.append(download_file[0])
                                full_path_file_list.append(full_path_file[0])
                                #print(''+self.__himawari_FTP+'/'+YYYY+MM+'/'+CHN+'/'+download_file[0]+'')
###    4km resolution band data
                            if down_4km:
                                geo_path = ''+self.__himawari_FTP+'/'+YYYY+MM+'/4KM/'+YYYY+''+MM+''+DD+''
                                for band_var in band4km:
                                    if (CHN == 'EXT' or CHN == 'VIS' or CHN == 'SIR') and (band_var == 'rad' or band_var == 'rfc' or band_var == 'rfy'):
                                        file_name = [''+YYYY+''+MM+''+DD+''+HH+''+MN+'.'+CHN.lower()+'.'+num+'.'+band_var+'.fld.4km.bin']
                                        download_file = [file_name[0] + '.bz2']
                                        full_path_file = [''+geo_path+'/'+download_file[0]+'']
                                        file_list.append(file_name[0])
                                        zip_file_list.append(download_file[0])
                                        full_path_file_list.append(full_path_file[0])
                                        #print(full_path_file[0])
                                    if CHN == 'TIR' and (band_var == 'rad' or band_var == 'tbb'):
                                        file_name = [''+YYYY+''+MM+''+DD+''+HH+''+MN+'.'+CHN.lower()+'.'+num+'.'+band_var+'.fld.4km.bin']
                                        download_file = [file_name[0] + '.bz2']
                                        full_path_file = [''+geo_path+'/'+download_file[0]+'']
                                        file_list.append(file_name[0])
                                        zip_file_list.append(download_file[0])
                                        full_path_file_list.append(full_path_file[0])
                                        #print(full_path_file[0])
###    4km resolution geo-info data for RGB images adjustment
            if geo_info:    
                        for geo_name in geo:
                            geo_path = ''+self.__himawari_FTP+'/'+YYYY+MM+'/4KM/'+YYYY+''+MM+''+DD+''
                            if geo_name == 'cap.flg':
                                file_name = [''+YYYY+''+MM+''+DD+''+HH+''+MN+'.'+geo_name+'.fld.bin']
                            else:
                                file_name = [''+YYYY+''+MM+''+DD+''+HH+''+MN+'.'+geo_name+'.fld.4km.bin']
                            download_file = [file_name[0] + '.bz2']
                            full_path_file = [''+geo_path+'/'+download_file[0]+'']
                            #print(full_path_file[0])
                            file_list.append(file_name[0])
                            zip_file_list.append(download_file[0])
                            full_path_file_list.append(full_path_file[0])
        return(file_list, zip_file_list, full_path_file_list) 

    def check_exist_sub_domain_file(self,full_path_file_list,remove_list_flag=True):
        import glob
        import netCDF4 as nc
        full_path_file_list = self.check_list(full_path_file_list)
        print('File numbers from generated list:',len(full_path_file_list))
        print('Check pre-processed sub-domain nc files...')
        for file_name in full_path_file_list:
            split_name = file_name.split('/')
            zip_file = split_name[-1]
            band_date, band, band_num  = self.name_info(file_name,convert2tbb=True)
            path_year = band_date[0:4]
            path_mon = band_date[4:6]
            nc_data_path = self.data_path + '/sub_domain_data/'+path_year+'/'+path_mon+''

            array_shape = self.band_array(band)
            nc_info = self.nc_file_info(file_name)
            lon_idx,lat_idx,local_lon,local_lat = self.lonlat_index(self.lon_range, self.lat_range, array_shape)
            unchecked_file = sorted(glob.glob(''+nc_data_path+'/'+nc_info[0]+'.nc'))
            if len(unchecked_file)>0 and remove_list_flag:
                nc_re = nc.Dataset(unchecked_file[0], 'r',  format='NETCDF4_CLASSIC')
                nclon_s = nc_re.variables['lon'][0]
                nclon_e = nc_re.variables['lon'][-1]
                nclat_s = nc_re.variables['lat'][0]
                nclat_e = nc_re.variables['lat'][-1]
                if local_lon[0] >= nclon_s and local_lon[-1] <= nclon_e:
                    if local_lat[0] >= nclat_s and local_lat[-1] <= nclat_e:
                        print('Data:',zip_file)
                        print('Already pre-processed to nc file')
                        print('Remove file name from file list')
                        full_path_file_list = list(filter(lambda x: file_name not in x, full_path_file_list)) 
        print('File numbers for downloading:',len(full_path_file_list))
        return(full_path_file_list) 

    def download(self,full_path_file_list,download_flag=True):
        import os
        import glob
        import wget
###    generate list if the input is string
        file_path = self.check_list(full_path_file_list)
        print('Downloading...')
        for download_file in file_path:
###    collect file information
            zip_file, file_name = self.name_info(download_file)        
            band_date, band, band_num  = self.name_info(file_name,convert2tbb=True)
            path_year = band_date[0:4]
            path_mon = band_date[4:6]
            download_path = self.data_path + '/compressed_data/'+path_year+'/'+path_mon+''
            os.makedirs(download_path, exist_ok=True)
###    check file already downloaded or not
            downloaded_file = sorted(glob.glob(download_path+'/'+zip_file)) 
            file_len = len(downloaded_file)
###    download file
            if download_flag and file_len < 1:
                print(''+zip_file+'')
                save_path = os.path.join(download_path, os.path.basename(download_file))
                try:
                    wget.download(''+download_file+'',save_path)
                except:
                   print('Data download failed')
            else: 
                print(''+zip_file+'','Data already exist')
###    check file downloaded or not
            downloaded_file = sorted(glob.glob(download_path+'/'+zip_file))
            file_len = len(downloaded_file)
            if file_len < 1:
                print(''+zip_file+'','no Data download')
                with open('no_file.txt', 'a') as file:
                    file.write(''+zip_file+'\n')
                file_path = list(filter(lambda x: file_name not in x, file_path))
        return(file_path)

    def unzip(self,file_name):
        import bz2
        import glob
###    generate list if the input is string
        file_name = self.check_list(file_name)        
        data_array = []
        output_file_list = []
        print('Extracting file')
        for download_file in file_name:
###    collect file information
            zip_file, output_file_name = self.name_info(download_file)
            band_date, band, band_num  = self.name_info(output_file_name,convert2tbb=True)
            path_year = band_date[0:4]
            path_mon = band_date[4:6]
            download_path = self.data_path + '/compressed_data/'+path_year+'/'+path_mon+''
###    check file downloaded or not 
            downloaded_file = sorted(glob.glob(download_path+'/'+zip_file))
            file_len = len(downloaded_file)
            if file_len > 0.5:
###    unzip file
                print(zip_file)
                bz2_file = bz2.BZ2File(download_path+'/'+zip_file,'rb')
                data = bz2_file.read()
                data_array.append(data)
                output_file_list.append(output_file_name)
        return(output_file_list,data_array)

    def generate_binary(self,file_name,data_array):
###    generate list if the input is string (or data_array is not a list)
        file_name = self.check_list(file_name)
        data_array = self.check_data_list(data_array)
        data_path = self.data_path
        print('Generating data')
        num = 0
        for output_file_name in file_name:
            print(output_file_name)
            with open(''+data_path+'/'+output_file_name+'','wb') as output:
                output.write(data_array[num])
                output.close()             
            num = num + 1

    def read_binary(self,file_name,data_array=[],binary_type=[]):
        import numpy as np
###    generate list if the input is string (or data_array is not a list)
        file_name = self.check_list(file_name)
        data_array = self.check_data_list(data_array)
        np_binary_data = []
        print('Read binary data in correct data type & array shape')
###    if there is input data_array 
        if len(data_array)>0:
            num = 0
            print('Loading data from memory')
            for output_file_name in file_name:
                #print(output_file_name)
                d_type, array_shape = self.binary_data_info(output_file_name, binary_type)
                binary_data = np.frombuffer(data_array[num],dtype=d_type).reshape(array_shape,array_shape)
                np_binary_data.append(binary_data)
                num = num + 1
            return(file_name,np_binary_data)
###    if there is no input data_array
        else:
            num = 0
            print('Loading data from binary')
            for output_file_name in file_name:
                #print(output_file_name)
                d_type, array_shape = self.binary_data_info(output_file_name, binary_type)
                binary_data = np.fromfile(''+self.data_path+'/'+output_file_name+'',dtype=d_type).reshape(array_shape,array_shape)
                np_binary_data.append(binary_data)
                num = num + 1
            return(file_name, np_binary_data)

    def convert(self,file_name,np_binary_array):
        import numpy as np
        file_name = self.check_list(file_name)
        np_binary_array = self.check_data_list(np_binary_array)
        output_band_tbb = []
        #print('Converting digits into Albedo or TBB')
        num = 0
        for output_file_name in file_name:
            print(output_file_name)
###    generate band date, band, band number
            band_date, band, band_num = self.name_info(output_file_name, convert2tbb=True)
            if band == '4km' or band == 'cap':
                print('4km datasets, physical variables already converted')
                output_band_tbb.append(np_binary_array[num])
            else:
###    LUT for converting count to tbb(albedo)
                print('Converting digits into Albedo or TBB')
                band_date = int(band_date)
                if band_date > 202212130449:
                    LUT_file = ['himawari_preparation/'+band+'.'+band_num+'.H09']
                elif band_date > 201802130249 and band_date < 201802140710:
                    LUT_file = ['himawari_preparation/'+band+'.'+band_num+'.H09']
                else:
                    LUT_file = ['himawari_preparation/'+band+'.'+band_num+'.H08']
            #print(LUT_file[0])
                LUT = np.loadtxt(LUT_file[0])
###    convert main process
                data_array = np_binary_array[num]
                valid_indices = (data_array >= 0) & (data_array < LUT.shape[0])
###    Initialize band_tbb with a default fill value
                band_tbb = np.full(data_array.shape, -999.0)
###    Map valid indices to Albedo or TBB values
                band_tbb[valid_indices] = LUT[data_array[valid_indices], 1]
                band_tbb = np.array(band_tbb,dtype='<f4')
                output_band_tbb.append(band_tbb)
            num = num + 1
        return(output_band_tbb)

    def sub_domain_extract(self,output_band_tbb):
        output_band_tbb = self.check_data_list(output_band_tbb)
        sub_domain_band_tbb = []
        sub_domain_lon = []
        sub_domain_lat = []
        print('Extracting sub-domain data')
        for i in range(0,len(output_band_tbb)):
            band_data = output_band_tbb[i]
            array_shape = band_data.shape
            #scale_factor = int(24000/array_shape[0])
            lon_idx,lat_idx,local_lon,local_lat = self.lonlat_index(self.lon_range,self.lat_range,array_shape[0])
            r_band_data = band_data[::-1]
            output_band_data = r_band_data[lat_idx[0]:lat_idx[1],lon_idx[0]:lon_idx[1]]
            sub_domain_band_tbb.append(output_band_data)
            sub_domain_lon.append(local_lon)
            sub_domain_lat.append(local_lat)
        return(sub_domain_band_tbb,sub_domain_lon,sub_domain_lat)

    def generate_nc(self,file_name,sub_domain_band_tbb,sub_domain_lon,sub_domain_lat,nc_output=True):
        import os
        import pickle
        import numpy as np
        from netCDF4 import Dataset
        from numpy import dtype

        print('Generating nc file')
        num = 0
        for output_file_name in file_name:
            output_name, band_date, var_name, long_name, units, missing_value = self.nc_file_info(output_file_name)
            path_year = band_date[0:4]
            path_mon = band_date[4:6]
            nc_data_path = self.data_path + '/sub_domain_data/'+path_year+'/'+path_mon+''
            os.makedirs(nc_data_path, exist_ok=True)
    
            band_data = sub_domain_band_tbb[num]
            local_lon = sub_domain_lon[num]
            local_lat = sub_domain_lat[num]
            array_size = band_data.shape
            if nc_output:
                tin = np.double(np.arange(0,2,1))
                ncout = Dataset(''+nc_data_path+'/'+output_name+'.nc', 'w', format='NETCDF4')
                ncout.createDimension('time', None)  # unlimited
                ncout.createDimension('lat', array_size[0])
                ncout.createDimension('lon', array_size[1])

               # create time axis
                time = ncout.createVariable('time', dtype('double').char, ('time',))
                time.long_name = 'UTC time'
                time.units = 'days since '+band_date[0:4]+'-'+band_date[4:6]+'-'+band_date[6:8]+' '+band_date[8:10]+':'+band_date[10:12]+':00'
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
                dout = ncout.createVariable(var_name, dtype('float').char, ('time', 'lat', 'lon'))
                dout.long_name = long_name
                dout.units = units
                dout.missing_value = missing_value
               # copy axis from original dataset
                time[:] = tin[0]
                lon[:] = local_lon[:]
                lat[:] = local_lat[:]
                output = np.zeros((1,array_size[0],array_size[1]))
                output[0,:,:] = band_data
                dout[:] = output[:]
                ncout.close()
                print(''+output_name+'.nc generated')
            else:
               # output .pkl file
                with open(''+output_name+'.pkl', 'wb') as f:
                    pickle.dump(band_data, f)
                print(''+output_name+'.pkl generated')
            num = num + 1 

    def name_info(self, file_name, read_binary=False, convert2tbb=False):
###    extract information from file list
        split_name = file_name.split('/')
###
        binary_types = self.string_info(binary_info=True)
        zip_file = split_name[-1]
        split_parts = zip_file.split('.')
        for b_type in binary_types:
            if b_type in split_parts:
                index = split_parts.index(b_type)
                unzip_file_name = ".".join(split_parts[:index + 1])
                binary_type = b_type
                break
###     
        band_types = self.string_info(band_info=True)
        for b_type in band_types:
            if b_type in split_parts:
                index = split_parts.index(b_type)
                band = b_type
                break
###
        band_num_types = self.string_info(band_num_info=True)
        for b_type in band_num_types:
            if b_type in split_parts:
                index = split_parts.index(b_type)
                band_num = b_type
                break
        #split_type = unzip_file_name.split('.')
        band_date = split_parts[0]
        if read_binary:
            return(binary_type, band)
        elif convert2tbb:
            return(band_date, band, band_num)
        else:
            return(zip_file, unzip_file_name)
    
    def band_array(self,band):
        if band == 'tir':
            array_shape = 6000
        elif band == 'sir':
            array_shape = 6000
        elif band == 'vis':
            array_shape = 12000
        elif band == 'ext':
            array_shape = 24000
        elif band == '4km' or band == 'cap':
            array_shape = 3000
        return(array_shape)

    def binary_data_info(self, output_file_name, binary_type):
        if len(binary_type) < 1:
            binary_type, band = self.name_info(output_file_name,read_binary=True)
            #print(binary_type)
        if binary_type == 'geoss':
            dtype  = '>u2'
            #print(dtype)
            array_shape = self.band_array(band)
        elif binary_type == 'bin':
            bin2array = '4km'
            dtype  = '>f4'
            if band == 'cap':
                dtype  = '>u2'
            array_shape = self.band_array(bin2array)
            #print(dtype)
        elif binary_type == 'dat':
            dtype  = '<f4'
            #print(dtype)
            array_shape = self.band_array(band)
        return(dtype,array_shape)
 
    def check_list(self,unchecked_var):
        if isinstance(unchecked_var, str):
            unchecked_var = [unchecked_var]
        return(unchecked_var)

    def check_data_list(self,unchecked_array):
        if isinstance(unchecked_array, list):
            unchecked_array = unchecked_array
        else:
            unchecked_array = [unchecked_array]
        return(unchecked_array)

    def string_info(self, binary_info=False, band_info=False, band_num_info=False,nc_info=False,nc_4km_info=False,nc_4km_var_info=False,time_info=False):
###
        binary_types = ['geoss', 'bin', 'dat']
        band_types = ['4km','cap', 'ext', 'vis', 'sir', 'tir']
        band_num_types = ['4km','cap', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10']
###    for nc file info
        band_table = ['ext01','vis01','vis02','vis03','sir01','sir02','tir01','tir02',
                      'tir03','tir04','tir05','tir06','tir07','tir08','tir09','tir10']
        AHI_band_num_table = ['03','01','02','04','05','06','13','14',
                              '15','16','07','08','09','10','11','12']
        data_type_table_4km = ['sun','sat','lat','lng','grd','cap','rad','rfc','rfy','tbb']
        angle_types = ['azm','zth']
###    for 4km nc var info
        nc_data_type_table_4km = ['sun_azm','sun_zth','sat_azm','sat_zth','lat','lng',
                                  'grd_time_mjd_hms','cap_flg','rad','rfc','rfy','tbb']
        var_name_table = ['sun_azm','sun_zth','sat_azm','sat_zth','latt','lng',
                          'grd_time_mjd_hms','cap_flg','rad','rfc','rfy','tbb']
        long_name_table = ['Solar azimuth angle (South direction is zero, clockwise rotation)', 
                           'Solar zenith angle', 
                           'Sensor azimuth angle (South direction is zero, clockwise rotation)', 
                           'Sensor zenith angle','Latitude','Longitude', 
                           'Scanning time（Normalized 0 to 1, i.e., 12:00 UTC is 0.5)', 
                           'Cloud flag (Daytime and over ocean only. More than 1 represents cloud)', 
                           'irradiance','spectral reflectance', 
                           'spectral reflectance','brightness temperature']
        units_table = ['degree','degree','degree','degree','degree','degree',
                       'None','None','W m-2 sr-1 μm-1','None','%','K']
        missing_value_table = [-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,
                               -99.0,-99.0,-99.0,-99.0,-99.0,-99.0]
        time_types = ['days','hours','minutes']
        if binary_info:
            return(binary_types)
        if band_info:
            return(band_types)
        if band_num_info:
            return(band_num_types)
        if nc_info:
            return(band_table, AHI_band_num_table)
        if nc_4km_info:
            return(band_table, AHI_band_num_table, data_type_table_4km, angle_types)
        if nc_4km_var_info:
            return(nc_data_type_table_4km,var_name_table,long_name_table,units_table,missing_value_table)
        if time_info:
            return(time_types)

    def nc_file_info(self, output_file_name):
###    name rule of nc file (band data & 4km data)
        band_date, band, band_num = self.name_info(output_file_name, convert2tbb=True)
        if band == 'cap':
            band = '4km'
        if band != '4km':
            band_table, AHI_band_num_table = self.string_info(nc_info=True)
            pos_idx = band_table.index(''+band+''+band_num+'')
            output_name = ''+ band_date +'_band_'+ AHI_band_num_table[pos_idx] +''
            var_name = 'tbb'
            long_name = 'Brightness Temperature(albedo)' 
            units = 'K or relf'
            missing_value = -999.
        else:
            band_table, AHI_band_num_table,data_type_table_4km,angle_types = self.string_info(nc_4km_info=True)
            nc_data_type_table_4km,var_name_table,long_name_table,units_table,missing_value_table = self.string_info(nc_4km_var_info=True)
            split_name = output_file_name.split('/')
            zip_file = split_name[-1]
            split_parts = zip_file.split('.')
###    generate 4km var name
            for b_type in data_type_table_4km:
                if b_type in split_parts:
                    index = split_parts.index(b_type)
                    data_type = b_type
                    if data_type in ('sun', 'sat'):
                        for a_type in angle_types: 
                            if a_type in split_parts:
                                data_type = ''+data_type+'_'+a_type+''
                                break
                    if data_type == 'grd':
                        data_type = 'grd_time_mjd_hms'
                    if data_type == 'cap':
                        data_type = 'cap_flg'
                    break
###    band data file name
            if data_type in ('rad', 'rfc', 'rfy', 'tbb'):
                band_types = self.string_info(band_info=True)
                band_types.pop(0)
                for b_type in band_types:
                    if b_type in split_parts:
                        index = split_parts.index(b_type)
                        band_4km = b_type
                        break
                band_num_types = self.string_info(band_num_info=True)
                band_num_types.pop(0)
                for b_type in band_num_types:
                    if b_type in split_parts:
                        index = split_parts.index(b_type)
                        band_num_4km = b_type
                        break
                pos_idx = band_table.index(''+band_4km+''+band_num_4km+'')
                output_name = ''+ band_date +'_4km_band_'+ AHI_band_num_table[pos_idx] +'_'+ data_type +''
###    geo data file name
            else:
                output_name = ''+ band_date +'_4km_'+ data_type +''
###    generate nc var info
            if data_type in nc_data_type_table_4km:
                pos_idx = nc_data_type_table_4km.index(data_type)
                var_name = var_name_table[pos_idx]
                long_name = long_name_table[pos_idx]
                units = units_table[pos_idx]
                missing_value = missing_value_table[pos_idx]
        return(output_name, band_date, var_name, long_name, units, missing_value)

    def lonlat_index(self, lon_range, lat_range, array_shape):
        import numpy as np
        x = np.arange(850024,2050024,50)/10000
        y = np.arange(-599976,600024,50)/10000
##
        ta = lon_range[0]
        diff = np.abs(x-ta)
        index = np.argmin(diff)
        x_dis_bot = np.mod(index,8)
###    lon array start
        lon_s = index - x_dis_bot

        ta = lon_range[1]
        diff = np.abs(x-ta)
        index = np.argmin(diff)
        x_dis_top = np.mod(index,8)
###    lon array end
        if x_dis_top > 0.5:
            lon_e = index + (8 - x_dis_top)
        else:
            lon_e = index
##
        ta = lat_range[0]
        diff = np.abs(y-ta)
        index = np.argmin(diff)
        y_dis_bot = np.mod(index,8)
###    lat array start
        lat_s = index - y_dis_bot

        ta = lat_range[1]
        diff = np.abs(y-ta)
        index = np.argmin(diff)
        y_dis_top = np.mod(index,8)
###    lat array end
        if y_dis_top > 0.5:
            lat_e = index + (8 - y_dis_top)
        else:
            lat_e = index
        scale_factor = int(24000/array_shape)
        resol = 50*scale_factor
        center = resol/2
        xx = np.arange(850000+center,2050000+center,resol)/10000
        yy = np.arange(-600000+center,600000+center,resol)/10000
#        print('lon',xx[0],xx[-1])
#        print('lat',yy[0],yy[-1])

        lon_idx = [int(lon_s/scale_factor),int(lon_e/scale_factor)]
        lat_idx = [int(lat_s/scale_factor),int(lat_e/scale_factor)]
#  print(lon_idx[0],lon_idx[1])
#  print(lat_idx[0],lat_idx[1])
        local_lon=xx[lon_idx[0]:lon_idx[1]]
        local_lat=yy[lat_idx[0]:lat_idx[1]]
#        print(local_lon[0],local_lon[-1])
#        print(local_lat[0],local_lat[-1])
        return(lon_idx,lat_idx,local_lon,local_lat)

    def information(self, detail=False):
        self.detail = detail
        print('Himawari-8 data strat from 2015/07/07 0200UTC.')
        print('Full disk covered area: 85E – 205E (155W), 60N – 60S')
        print('')
        print('Band data & geo-info data that can be downloaded')
        print('-------------------------------------------------------')
        print('[EXT] 01:Band03')
        print('[VIS] 01:Band01 02:Band02 03:Band04')
        print('[SIR] 01:Band05 02:Band06')
        print('[TIR] 01:Band13 02:Band14 03:Band15 04:Band16 05:Band07')
        print('      06:Band08 07:Band09 08:Band10 09:Band11 10:Band12')
        print('[GEO] Solar azimuth angle(sun.azm) Solar zenith angle(sun.zth)')
        print('      Sensor azimuth angle(sat.azm) Sensor zenith angle(sat.zth)')
        print('-------------------------------------------------------')
        print('')
        print('Download period & data template')
        print('-------------------------------------------------------')
        print('year = [\'2015\']               # Year:   from 2015')
        print('mon  = [\'07\',\'09\']            # Month:  01 02 ... 12')
        print('day  = [\'07\']                 # Day:    01 02 ... 30 31')
        print('hour = [\'02\']                 # Hour:   00 01 ... 23')
        print('minn = [\'00\']                 # Minute: 00 10 20 30 40 50')
        print('Set band type')
        print('band = ['+'EXT'+', '+'VIS'+', '+'SIR'+', '+'TIR'+']   # band: VIS TIR SIR EXT')
        print('band_num = ['+'1,2,3'+']            # band number: 1 2 ... 10')
        print('geo = [\'sun.azm\', \'sun.zth\']  # geo-info: sun.azm sun.zth sat.azm sat.zth ...')
        print('-------------------------------------------------------')
        print('You don\'t need to download a continuous time interval.')
        print('Instead, you can choose a customized time resolution based on your research purpose:')
        print('Such as one entry per hour, one entry per day, ... or even one entry per year.')
        print('-------------------------------------------------------')
        print('')
        print('[EXT], [VIS], [SIR]: albedo [%]; [TIR]: brightness temperature [K]')
        print('')
        print('[EXT] Resolution: 0.005 degree;   array_size: (24000,24000)')
        print('[VIS] Resolution: 0.01 degree;    array_size: (12000,12000)')
        print('[SIR] Resolution: 0.02 degree;    array_size: (6000,6000)')
        print('[TIR] Resolution: 0.02 degree;    array_size: (6000,6000)')
        print('[GEO] Resolution: 0.04 degree;    array_size: (3000,3000)')
        print('')
        print('Reference: http://www.cr.chiba-u.jp/databases/GEO/H8_9/FD/index.html')
        if detail:
            print('')
            print('Extra information')
            print('-------------------------------------------------------')
            print('GrADS control file lat lon format')
            print('[EXT]')
            print('xdef 24000 linear 85.0025  0.005')
            print('ydef 24000 linear -59.9975 0.005')
            print('[VIS]')
            print('xdef 12000 linear 85.005  0.01')
            print('ydef 12000 linear -59.995 0.01')
            print('[SIR] & [TIR]')
            print('xdef 6000 linear 85.01  0.02')
            print('ydef 6000 linear -59.99 0.02')
            print('[GEO]')
            print('xdef 3000 linear 85.02  0.04')
            print('ydef 3000 linear -59.98 0.04')
            print('')
            print('-------------------------------------------------------')
            print('Full geometries dataset include:')
            print('[GEO]')
            print('Solar azimuth angle(sun.azm) Solar zenith angle(sun.zth)')
            print('Sensor azimuth angle(sat.azm) Sensor zenith angle(sat.zth)')
            print('Latitude(lat) Longitude(lng)')
            print('Scanning time(grd.time.mjd.hms) Cloud flag(cap)')
            print('')
            print('All band data also provide 0.04degree(4km) resolution physical variables converted data:')
            print('[BAND]')
            print('YYYYMMDDHHMN.xxx.ZZ.rad.fld.4km.bin.bz2 (xxx: ext, vis, sir, tir; ZZ: CEReS gridded data band number) ')
            print('ext, vis, sir, tir irradiance (unit: W m-2 sr-1 μm-1) ')
            print('YYYYMMDDHHMN.xxx.ZZ.rfc.fld.4km.bin.bz2 (xxx: ext, vis, sir; ZZ: CEReS gridded data band number) ')
            print('ext, vis, sir spectral reflectance (dimensionless) ')
            print('YYYYMMDDHHMN.xxx.ZZ.rfy.fld.4km.bin.bz2 (xxx: ext, vis, sir; ZZ: CEReS gridded data band number) ')
            print('ext, vis, sir spectral reflectance (%)')
            print('YYYYMMDDHHMN.tir.ZZ.tbb.fld.4km.bin.bz2 (ZZ: CEReS gridded data band number) ')
            print('tir (only) brightness temperature (Tbb (K)) ')
            print('-------------------------------------------------------')

    def generate_time_list(self,time_period=[],time_delta=[],year=[],mon=[],day=[],hour=[],minn=[]):
        self.year = year
        self.mon  = mon
        self.day  = day
        self.hour = hour
        self.minn = minn
      
###  two types for generating the time list
###  one for continious time range; the other for specific time
        time_list = [] 
###  continious time range
        if len(time_period) > 0:
            self.time_period = time_period
            from datetime import datetime, timedelta
            time_types = self.string_info(time_info=True)
            date_1 = time_period[0]
            date_2 = time_period[1]
            s_date = datetime(int(date_1[0:4]),int(date_1[4:6]),int(date_1[6:8]))
            e_date = datetime(int(date_2[0:4]),int(date_2[4:6]),int(date_2[6:8]))
            delta_value = [0,0,0]
            num = 0
            for t_type in time_types:
                for time_part in time_delta:
                    split_parts = time_part.split('=')
                    if t_type in split_parts:
                        delta_value[num] = int(split_parts[1])
                num = num + 1
            delta = timedelta(days=delta_value[0],hours=delta_value[1],minutes=delta_value[2])
            time_range = []
            if sum(delta_value)>0:
                while s_date <= e_date:
                    time_range.append(s_date)
                    s_date += delta
            else:
                time_range.append(s_date)

            for date in time_range:
                str_date = str(date)
                time = str_date[0:4] + str_date[5:7] + str_date[8:10] + str_date[11:13] + str_date[14:16]
                time_list.append(time)

###  specific time
        else:
          for YYYY in year:
            for MM in mon:
              for DD in day:
                for HH in hour:
                  for MN in minn:
                    time = YYYY + MM + DD + HH + MN
                    time_list.append(time) 
###
        return(time_list)
