from satellite_preprocess import Himawari

d_path = '/data/C.jerryjerry9/hima_download/himawari_data'
path = '/data/C.jerryjerry9/hima_download'
lat = [-10, 50]
lon = [90, 180]
time_period = ['20190123','20190124']


himawari = Himawari.Himawari(work_path=path,data_path=d_path,lat_range=lat,lon_range=lon,time_period=time_period)
#print(himawari.work_path)
#print(himawari.lat_range)
#print(himawari.lon_range)
#print(himawari.time_period)
#print(himawari.data_path)

himawari.information()
