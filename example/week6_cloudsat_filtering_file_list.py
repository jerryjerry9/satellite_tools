from ael_satellite_tools.preprocess import CloudSat

lat = [-10, 60]
lon = [90, 150]

cloudsat = CloudSat(work_path=[],lat_range=lat,lon_range=lon)

time_period = ['2016122000','2016122624']

product_name = 'hdf-GEOPROF'
full_path_file_list_geo = cloudsat.generate_list(product_name, time_period)

product_name = 'hdf-GEOPROF-LIDAR'
full_path_file_list_lidar = cloudsat.generate_list(product_name, time_period)

match_file_list_geo,match_file_list_lidar = \
cloudsat.cross_product_match(full_path_file_list_geo,
                             full_path_file_list_lidar)

sub_domain_geo = cloudsat.sub_domain_check(match_file_list_geo,mask_output=False)
sub_domain_lidar = cloudsat.sub_domain_check(match_file_list_lidar,mask_output=False)

label = True
cloudsat.plot_track(sub_domain_geo,coast_line_color='olive',
                    lonlat_step=4,font_size=24,loc='left',pad=40,
                    figure_title='CloudSat Track Demo', 
                    prefix='cloudsat_track_demo', 
                    figure_path='cloudsat_fig', dpi=300,
                    granule_label=False,id_label=label,
                    date_label=label,utc_label=label,save_fig=True)

