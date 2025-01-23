###
# Satellite process tools
# it should include: 
# download; pre-process; plotting;
# pre-process may include extract sub-domain, re-grid, composite ..... 
###


class Preparation:
    def __init__(self,work_path,lat_range,lon_range,time_period):
        self.work_path = work_path
        self.lat_range = lat_range
        self.lon_range = lon_range
        self.time_period = time_period

class Himawari(Preparation):
    def __init__(self,work_path=None,lat_range=[-10,50],lon_range=[90,180],time_period=['20150707'],data_path='/data/dadm1/obs/Himawari'):
        super().__init__(work_path,lat_range,lon_range,time_period)
        self.data_path = data_path

    def test_part(self):
        data_path = self.data_path
        print(data_path)
        print(self.data_path)

    def generated_list(self,time_period=[],year=[],mon=[],day=[],hour=[],minn=[],band=[],band_num=[],band4km=[],geo=[]):
###  prepare vars
        self.__himawari_FTP = 'ftp://hmwr829gr.cr.chiba-u.ac.jp/gridded/FD/V20190123'
        self.time_period = time_period
        self.year = year
        self.mon  = mon
        self.day  = day
        self.hour = hour
        self.minn = minn
        self.band = band
        self.band_num = band_num
        self.geo  = geo
        self.band4km = band4km
###  generating dwonload flag
        if len(band) > 0 and len(band_num) > 0:      
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
        if len(time_period) > 0:
           continuous_range = True
           year_list, mon_list , day_list = self.generate_day_range(time_period=time_period,continuous_range=continuous_range)
        else:
           continuous_range = False           
           year_list, mon_list , day_list = self.generate_day_range(year=year,mon=mon,day=day,continuous_range=continuous_range)
###  generating file list for download from FTP           
        file_list = []
        zip_file_list = [] 
        full_path_file_list = []
###    time period
        for time_num in range(0,len(year_list)):
          YYYY = year_list[time_num]
          MM = mon_list[time_num]
          DD = day_list[time_num]
          for HH in hour:
            for MN in minn:
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
                  file_name = [''+YYYY+''+MM+''+DD+''+HH+''+MN+'.'+geo_name+'.fld.4km.bin']
                  download_file = [file_name[0] + '.bz2']
                  full_path_file = [''+geo_path+'/'+download_file[0]+'']
                  #print(full_path_file[0])
                  file_list.append(file_name[0])
                  zip_file_list.append(download_file[0])
                  full_path_file_list.append(full_path_file[0])
        return(file_list, zip_file_list, full_path_file_list) 



    def download(self,file_path):
        import os
        import glob
        import wget
###    generate list if the input is string
        file_path = self.check_list(file_path)
        print('Downloading...')
        for download_file in file_path:
###    collect file information
            zip_file, file_name = self.name_info(download_file)
###    download file        
            print(''+zip_file+'') 
            #wget.download(''+download_file+'')
###    check file downloaded or not
            downloaded_file = sorted(glob.glob(zip_file))
            file_len = len(downloaded_file)
            if file_len < 1:
                with open('no_file.txt', 'a') as file:
                    file.write(''+zip_file+'\n')



    def unzip(self,file_name):
        import bz2
        import glob
        print('Extract file')
###    generate list if the input is string
        file_name = self.check_list(file_name)        
        data_array = []
        output_file_list = []
        for download_file in file_name:
###    collect file information
            zip_file, output_file_name = self.name_info(download_file)
###    check file downloaded or not 
            downloaded_file = sorted(glob.glob(zip_file))
            file_len = len(downloaded_file)
            if file_len > 0.5:
###    unzip file
                print(zip_file)
                bz2_file = bz2.BZ2File(zip_file,'rb')
                data = bz2_file.read()
                data_array.append(data)
                output_file_list.append(output_file_name)
        return(output_file_list,data_array)



    def output_binary(self,file_name,data_array):
###    generate list if the input is string (or data_array is not a list)
        file_name = self.check_list(file_name)
        data_array = self.check_data_list(data_array)
        print('Generating data')
        data_path = self.data_path
        num = 0
        for output_file_name in file_name:
            print(output_file_name)
            with open(''+data_path+'/'+output_file_name+'','wb') as output:
                output.write(data_array[num])
                output.close()             
            num = num + 1

    def read_binary(self):
        pass


    def convert(self,file_name,data_array,file_list=False):
        if len(data_array)>0:
            print('convert from memory')
        else:
            print('convert from file')

    def sub_domain(self):
        pass

    def output_nc(self):
        pass 

    def name_info(self,file_name):
###    extract information from file list
        split_name = file_name.split('/')
###
        binary_types = ['geoss', 'bin', 'dat']
        zip_file = split_name[-1]
        split_parts = zip_file.split('.')
        for b_type in binary_types:
            if b_type in split_parts:
                index = split_parts.index(b_type)
                unzip_file_name = ".".join(split_parts[:index + 1])
                break
###     
        split_type = unzip_file_name.split('.')
        binary_type = split_type[-1]
        return(zip_file,unzip_file_name)
    
    def band_array(self,band):
        if band == 'tir':
            array_shape = 6000
        elif band == 'sir':
            array_shape = 6000
        elif band == 'vis':
            array_shape = 12000
        elif band == 'ext':
            array_shape = 24000
        return(array_shape)
    
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


    def generate_nc_process(self,file_list,unzip_flag=True,convert_flag=True,output_flag=False,rm_flag=True):
###    unzip file
#                if unzip_flag:
#                    print('Extract file')
#                    print(file_name)
#                    if output_flag and data_type == 'bin':
#                        self.unzip(file_name,output_data=output_flag)
#                        print('Generate 4km .bin file')
#                    elif output_flag and data_type == 'geoss':
#                        data_array = self.unzip(file_name,output_data=output_flag)
#                        if convert_flag:
#                            data_array = self.unzip(file_name,output_data=False)
#                            print('Convert and generate .dat file')
#                        else:
#                            self.unzip(file_name,output_data=output_flag)
#                            print('Gerenate .geoss file ... Do you realy need this?')
#                    else:
#                        data_array = self.unzip(file_name,output_data=output_flag)
#                        print(len(data_array))
#                        if convert_flag and data_type == 'geoss':
#                            print('convert!')
###    remove .bz2 files
#                if rm_flag:
#                    os.remove(zip_file)

        pass


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

    def generate_day_range(self,time_period=[],year=[],mon=[],day=[],continuous_range=False):
###  two types for generating the day range list
###  one for continious day range; the other for specific day list
        year_list = []
        mon_list = []
        day_list = []
###  continious day range
        if continuous_range: 
          from datetime import datetime, timedelta
          self.time_period = time_period
          date_1 = time_period[0]
          date_2 = time_period[1]
          s_date = datetime(int(date_1[0:4]),int(date_1[4:6]),int(date_1[6:8]))
          e_date = datetime(int(date_2[0:4]),int(date_2[4:6]),int(date_2[6:8])) + timedelta(days=1)
          delta = timedelta(days=1)
          #delta = timedelta(hours=5)
          #delta = timedelta(minutes=10)
          time_range = []

          while s_date < e_date:
            time_range.append(s_date)
            s_date += delta

          for date in time_range:
            str_date = str(date)
            year_list.append(str_date[0:4])
            mon_list.append(str_date[5:7])
            day_list.append(str_date[8:10])
###  specific day list
        else:
          self.year = year
          self.mon = mon
          self.day = day
          for YYYY in year:
            for MM in mon:
              for DD in day:
                year_list.append(YYYY)
                mon_list.append(MM)
                day_list.append(DD)
###
        return(year_list, mon_list, day_list)

