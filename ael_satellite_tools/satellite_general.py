"""
satellite_general.py

This module provides general functions for satellite associated modules

"""

class Preparation:
    def __init__(self,work_path,lat_range,lon_range,time_period):
        self.work_path = work_path
        self.lat_range = lat_range
        self.lon_range = lon_range
        self.time_period = time_period

    def generate_time_list(self,time_period=[],time_delta=[],year=[],mon=[],day=[],hour=[],minn=[]):
        self.year = year
        self.mon  = mon
        self.day  = day
        self.hour = hour
        self.minn = minn

###  two types for generating the time list
###  one for continuous time range; the other for specific time
        time_list = []
###  continuous time range
        if len(time_period) > 0:
            self.time_period = time_period
            from datetime import datetime, timedelta
            time_types = self.string_info(time_info=True)
            date_1 = time_period[0]
            date_2 = time_period[1]
            s_date = datetime(int(date_1[0:4]),int(date_1[4:6]),int(date_1[6:8]),int(date_1[8:10]),int(date_1[10:12]))
            e_date = datetime(int(date_2[0:4]),int(date_2[4:6]),int(date_2[6:8]),int(date_1[8:10]),int(date_1[10:12]))
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

    def check_list(self,unchecked_var):
        if isinstance(unchecked_var, list):
            unchecked_var = unchecked_var
        else: 
            unchecked_var = [unchecked_var]
        return(unchecked_var)

    def check_data_list(self,unchecked_array):
        if isinstance(unchecked_array, list):
            unchecked_array = unchecked_array
        else:
            unchecked_array = [unchecked_array]
        return(unchecked_array)

    def calculate_julian_date(self, year, month, day):
        from datetime import datetime
        date = datetime(year, month, day)
        julian_date = date.toordinal() + 1721425.5
        return(julian_date)

    def read_nc_boundary(self,file_name,lon_range,lat_range,array_shape):
        import netCDF4 as nc
        local_lon,local_lat = self.lonlat_index(lon_range, lat_range, array_shape=array_shape,lonlat_only=True)
        nc_re = nc.Dataset(file_name[0], 'r',  format='NETCDF4_CLASSIC')
        nclon_s = nc_re.variables['lon'][0]
        nclon_e = nc_re.variables['lon'][-1]
        nclat_s = nc_re.variables['lat'][0]
        nclat_e = nc_re.variables['lat'][-1]
        ## for convert back to lon lat range
        nclon_2 = nc_re.variables['lon'][1]
        nclat_2 = nc_re.variables['lat'][1]
        step = nclon_2 - nclon_s
        nc_lon_range = [nclon_s]
        bot_lon = int((int(round(nclon_s*10000)) - int(round(step*10000)/2))/10000)
        top_lon = int((int(round(nclon_e*10000)) + int(round(step*10000)/2))/10000)
        bot_lat = int((int(round(nclat_s*10000)) - int(round(step*10000)/2))/10000)
        top_lat = int((int(round(nclat_e*10000)) + int(round(step*10000)/2))/10000)
        nc_lon_range = [bot_lon,top_lon]
        nc_lat_range = [bot_lat,top_lat]


        if local_lon[0] >= nclon_s and local_lon[-1] <= nclon_e and local_lat[0] >= nclat_s and local_lat[-1] <= nclat_e:
           domain_flag = True
        else:
           domain_flag = False
        return (domain_flag,nc_lon_range,nc_lat_range)

    def extend_lonlat(self,nc_lon_range,nc_lat_range):
        import numpy as np
        new_lon_range = []
        new_lat_range = []
        lat_min = np.minimum(nc_lat_range[0], self.lat_range[0])
        lon_min = np.minimum(nc_lon_range[0], self.lon_range[0])
        lat_max = np.maximum(nc_lat_range[1], self.lat_range[1])
        lon_max = np.maximum(nc_lon_range[1], self.lon_range[1])
        new_lon_range.append(lon_min) 
        new_lon_range.append(lon_max)      
        new_lat_range.append(lat_min)
        new_lat_range.append(lat_max)
        return(new_lon_range, new_lat_range)

    def lonlat_index(self, lon_range, lat_range, array_shape=[], resolution=[], lonlat_only=False):
        import numpy as np
        resolution = self.check_list(resolution)
        if len(resolution) > 0:
            array_shape = 12000/resolution[0]
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
        if lonlat_only:
            return(local_lon,local_lat)
        else:
            return(lon_idx,lat_idx,local_lon,local_lat)

    def band_name_convert(self,band_info,AHI2CEReS=True):
        band_info = self.check_list(band_info)
        band_table, AHI_band_num_table = self.string_info(nc_info=True)
        CEReS_type = []
        CEReS_num = []
        if AHI2CEReS:
            for band in band_info:
                AHI_num = str(band)
                if band < 10:
                   AHI_num = '0' + AHI_num
                pos_idx = AHI_band_num_table.index(AHI_num) 
                band_name = band_table[pos_idx]
                band_type = band_name[0:3].upper()
                band_num = int(band_name[3:5])
                if band_type not in CEReS_type:
                    CEReS_type.append(band_type)
                if band_num not in CEReS_num:
                    CEReS_num.append(band_num)
        return(CEReS_type,sorted(CEReS_num))

    def string_info(self, binary_info=False, band_info=False, band_num_info=False,nc_info=False,nc_4km_info=False,nc_4km_var_info=False,time_info=False,nc_plotting_list=False):
###
        binary_types = ['geoss', 'bin', 'dat']
        band_types = ['4km','cap', 'ext', 'vis', 'sir', 'tir']
        band_num_types = ['4km','cap', '01', '02', '03', '04', '05', '06', '07', '08', '09', '10']
###    for nc file info
        band_table = ['ext01','vis01','vis02','vis03','sir01','sir02','tir01','tir02',
                      'tir03','tir04','tir05','tir06','tir07','tir08','tir09','tir10']
        AHI_band_num_table = ['03','01','02','04','05','06','13','14',
                              '15','16','07','08','09','10','11','12']
        band_resolution = [24000, 12000, 12000, 12000, 6000, 6000, 6000, 6000,
                           6000, 6000, 6000, 6000, 6000, 6000, 6000, 6000]

        data_type_table_4km = ['sun','sat','lat','lng','grd','cap','rad','rfc','rfy','tbb']
        angle_types = ['azm','zth']
###    for 4km nc var info
        geo_name_table = ['sun.azm','sun.zth','sat.azm','sat.zth','lat','lng',
                                  'grd.time.mjd.hms','cap.flg']
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
        if nc_plotting_list:
            return(band_table, AHI_band_num_table, band_resolution, geo_name_table, nc_data_type_table_4km)


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

