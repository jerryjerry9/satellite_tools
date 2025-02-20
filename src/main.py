from satellite_preprocess import Himawari

## prepare enviromnent
path = []
data_path = '/data/C.jerryjerry9/satellite_tools/data/Himawari'
lat = [-10, 50]
lon = [90, 180]

himawari = Himawari.Himawari(work_path=path,data_path=data_path,lat_range=lat,lon_range=lon)
print(himawari.work_path)
print(himawari.lat_range)
print(himawari.lon_range)
print(himawari.time_period)
print(himawari.data_path)

## Set download time period
## two ways to set time period
## mehtod 1: self-defined specific download time
year = ['2017']   # Year: from 2015
mon = ['02']      # Month: 01 02 ... 12
day = ['10']      # Day: 01 02 ... 30 31
hour = ['00']     # Hour: 00 01 ... 23
minn = ['00']     # Minute: 00 10 20 30 40 50
## method 2: download data within a time range
time_period = ['201702100000','201702230000'] # YYYYMMDDHHMN ~ YYYYMMDDHHMN
time_delta = ['days=0','hours=1','minutes=0'] # time interval (minimum time resolution: 10 minutes)
## Set band type
band = ['TIR']  # band: VIS TIR SIR EXT
band_num = [1,2]  # band number: 1 2 ... 10 (digits)
geo = ['sun.azm','sun.zth','sat.azm','sat.zth'] # geo: sun.azm, sun.zth, sat.azm, sat.zth
band4km = [] # band_4km: rad, rfc, rfy, tbb

### downloading & pre-processing main body
## method 1: self-defined specific download time
file_list, zip_file_list, full_path_file_list = himawari.generate_list(year=year,mon=mon,day=day,hour=hour,minn=minn,band=band,band_num=band_num,geo=geo, band4km=band4km)

## method 2: download data within a time range
## Given time_period & time_delta will prioritize using method 2
#file_list, zip_file_list, full_path_file_list = himawari.generate_list(time_period=time_period,time_delta=time_delta,year=year,mon=mon,day=day,hour=hour,minn=minn,band=band,band_num=band_num,geo=geo, band4km=band4km)

himawari.pre_process(full_path_file_list,remove_list_flag=True)
#himawari.information()
