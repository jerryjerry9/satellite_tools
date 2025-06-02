"""
satellite_tools.py

This module provides functions for downloading & pre-processing satellite products
pre-processing may include extract sub-domain, re-grid, composite ...

planning satellite product: Himawari, CloudSat, RO, ...

"""
from ..satellite_general import Preparation


class Data_group:
    def __init__(self,band_data):
        self.band_data = band_data


class Himawari(Preparation):
    def __init__(self,work_path=None,lat_range=[-10,50],lon_range=[90,180],time_period=['20150707'],data_path='/data/dadm1/obs/Himawari'):
        from pathlib import Path 
        super().__init__(work_path,lat_range,lon_range,time_period)
        if self.work_path == None or self.work_path == []:
            self.work_path = Path().cwd()
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
        file_list, zip_file_list, ftp_path_file_list = self.generate_data_list(time_list=time_list,band=band,band_num=band_num,geo=geo,band4km=band4km,separate_ori_4km=separate_ori_4km)
        return (file_list, zip_file_list, ftp_path_file_list)

    def pre_process(self, ftp_path_file_list, remove_list_flag=True, download_flag=True):
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
        ftp_path_file_list = self.check_exist_sub_domain_file(ftp_path_file_list,remove_list_flag=remove_list_flag)

        for file_name in ftp_path_file_list:

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
        ftp_path_file_list = []
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
                                continue
                            elif CHN == 'SIR' and NUM > 2:
                                continue
                            elif CHN == 'EXT' and NUM > 1:
                                continue
                            num = str(NUM)
                            if NUM < 10:
                                num = '0'+ num
###    ori-resolution band data
                            if down_ori:
                                file_name = [''+YYYY+''+MM+''+DD+''+HH+''+MN+'.'+CHN.lower()+'.'+num+'.fld.geoss']
                                download_file = [file_name[0] + '.bz2']
                                full_path_file = [''+self.__himawari_FTP+'/'+YYYY+MM+'/'+CHN+'/'+download_file[0]+'']
                                file_list.append(file_name[0])
                                zip_file_list.append(download_file[0])
                                ftp_path_file_list.append(full_path_file[0])
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
                                        ftp_path_file_list.append(full_path_file[0])
                                        #print(full_path_file[0])
                                    if CHN == 'TIR' and (band_var == 'rad' or band_var == 'tbb'):
                                        file_name = [''+YYYY+''+MM+''+DD+''+HH+''+MN+'.'+CHN.lower()+'.'+num+'.'+band_var+'.fld.4km.bin']
                                        download_file = [file_name[0] + '.bz2']
                                        full_path_file = [''+geo_path+'/'+download_file[0]+'']
                                        file_list.append(file_name[0])
                                        zip_file_list.append(download_file[0])
                                        ftp_path_file_list.append(full_path_file[0])
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
                            ftp_path_file_list.append(full_path_file[0])
        return(file_list, zip_file_list, ftp_path_file_list) 

    def check_exist_sub_domain_file(self,ftp_path_file_list,remove_list_flag=True,extend_lonlat_range=True):
        import glob
        ftp_path_file_list = self.check_list(ftp_path_file_list)
        print('File numbers from generated list:',len(ftp_path_file_list))
        print('Check pre-processed sub-domain nc files...')
        if extend_lonlat_range == True:
            print('lon lat range will combine the nc file range & self-defined range')
        if extend_lonlat_range == False:
            print('lon lat range will use self-defined range')
        for file_name in ftp_path_file_list:
            split_name = file_name.split('/')
            zip_file = split_name[-1]
            band_date, band, band_num  = self.name_info(file_name,convert2tbb=True)
            path_year = band_date[0:4]
            path_mon = band_date[4:6]
            nc_data_path = self.data_path + '/sub_domain_data/'+path_year+'/'+path_mon+''

            array_shape = self.band_array(band)
            nc_info = self.nc_file_info(file_name)
            unchecked_file = sorted(glob.glob(''+nc_data_path+'/'+nc_info[0]+'.nc'))
            if len(unchecked_file)>0 and remove_list_flag:
                domain_flag,nc_lon_range,nc_lat_range = self.read_nc_boundary(unchecked_file,self.lon_range,self.lat_range, array_shape)
                if domain_flag == True:
                    print('Data:',zip_file)
                    print('Already pre-processed to nc file')
                    print('Remove file name from file list')
                    ftp_path_file_list = list(filter(lambda x: file_name not in x, ftp_path_file_list)) 
                else:
                    if extend_lonlat_range:
                        new_lon_range, new_lat_range = self.extend_lonlat(nc_lon_range,nc_lat_range)
                        self.lon_range = new_lon_range
                        self.lat_range = new_lat_range
        print('File numbers for downloading:',len(ftp_path_file_list))
        return(ftp_path_file_list) 

    def download(self,ftp_path_file_list,download_flag=True):
        import os
        import glob
        import wget
###    generate list if the input is string
        file_path = self.check_list(ftp_path_file_list)
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

    def read_binary(self,file_name,data_array=[],binary_type=[],band=[]):
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
                d_type, array_shape = self.binary_data_info(output_file_name, binary_type, band)
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
                d_type, array_shape = self.binary_data_info(output_file_name, binary_type, band)
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
                LUT = self.load_LUT(band_date,band,band_num)
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
            lon_idx,lat_idx,local_lon,local_lat = self.lonlat_index(self.lon_range,self.lat_range,array_shape=array_shape[0])
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

    def binary_data_info(self, output_file_name, binary_type, band):
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

    def load_LUT(self,band_date,band,band_num):
        from pathlib import Path
        import numpy as np
        current_dir = Path(__file__).parent
        lookup_path = str(current_dir / 'himawari_LUT')
        band_date = int(band_date)
        if band_date > 202212130449:
            LUT_file = [''+lookup_path+'/'+band+'.'+band_num+'.H09']
        elif band_date > 201802130249 and band_date < 201802140710:
            LUT_file = [''+lookup_path+'/'+band+'.'+band_num+'.H09']
        else:
            LUT_file = [''+lookup_path+'/'+band+'.'+band_num+'.H08']
        LUT = np.loadtxt(LUT_file[0])
        return(LUT)
