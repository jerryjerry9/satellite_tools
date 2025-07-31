"""
FS3 & FS7 RO profiles pre-process module
"""

from ..satellite_general import Preparation
from ..satellite_plotting import Plotting_tools


class RO(Preparation,Plotting_tools):
    def __init__(self,work_path=None,
                 lat_range=[-10,60],lon_range=[90,150],time_period=['200606'],
                 data_path='/data/dadm1/obs/RO_profile'):
        from pathlib import Path
        super().__init__(work_path,lat_range,lon_range,time_period)
        if self.work_path == None or self.work_path == []:
            self.work_path = Path().cwd()
        self.data_path = data_path
        self.date_range = []

    def generate_list(self,time_period,satellite_overlap='both'):
        import netCDF4 as nc
        import glob
        # 
        satellite_overlap = satellite_overlap.lower()
        fs3_file_path = ''+self.data_path+'/fs3'
        fs7_file_path = ''+self.data_path+'/fs7'
        # list = fs3 fs7 both
        time_period = self.check_list(time_period)
        year_list, ju_day_list, start_end_hour_list = self.convert_input_period(time_period)
        search_period = self.generate_period(year_list, ju_day_list, start_end_hour_list)
        full_path_file_list = []
        s_date = time_period[0][0:4]+'/'+time_period[0][4:6]+'/'+time_period[0][6:8]+' ' +time_period[0][8:10]+'UTC'
        e_date = time_period[1][0:4]+'/'+time_period[1][4:6]+'/'+time_period[1][6:8]+' ' +time_period[1][8:10]+'UTC'
        self.date_range = [s_date,e_date]
        for nc_time in search_period:
            filter_year = int(nc_time[0:4])
            filter_day =  int(nc_time[4:7])
            if filter_year == 2019:
                if filter_day <=344 and filter_day >= 274:
                   if satellite_overlap == 'both' or satellite_overlap == 'fs3':
                       search_files = self.ro_file_list(fs3_file_path,nc_time)
                       full_path_file_list.extend(search_files)
                   if satellite_overlap == 'both' or satellite_overlap == 'fs7':
                       search_files = self.ro_file_list(fs7_file_path,nc_time)
                       full_path_file_list.extend(search_files)
            else:
                search_files = self.ro_file_list(fs3_file_path,nc_time)
                full_path_file_list.extend(search_files)
                search_files = self.ro_file_list(fs7_file_path,nc_time)
                full_path_file_list.extend(search_files)
        return(full_path_file_list)
    
    def ro_file_list(self,file_path,nc_time):
###     wetPf2_C2E3.2024.144.21.30.G24_0001.0001_nc
        import glob
        year = nc_time[0:4]
        j_day =  nc_time[4:7]
        utc_time = nc_time[7:9]
        product_path = file_path +'/'+ year +'/'+ j_day + '/'
        if len(nc_time)>7:
            utc_time = nc_time[7:9]
            target_file = product_path + '*'+ ''+year+'.'+j_day+'.'+utc_time+'.*_nc'
        else:
            target_file = product_path + '*'+ ''+year+'.'+j_day+'.*_nc'
        search_files = sorted(glob.glob(target_file))
        return(search_files)

    def sub_domain_check(self,full_path_file_list,
                         extracted_lon_range=[],extracted_lat_range=[],
                         lonlat_list=False):
        import netCDF4 as nc
        import numpy as np
        full_path_file_list = self.check_list(full_path_file_list)
###  lon lat range for checking
        if len(extracted_lon_range)<1 and len(extracted_lat_range)<1:
            lon_s = self.lon_range[0]
            lon_e = self.lon_range[1]
            lat_s = self.lat_range[0]
            lat_e = self.lat_range[1]
        else:
            lon_s = extracted_lon_range[0]
            lon_e = extracted_lon_range[1]
            lat_s = extracted_lat_range[0]
            lat_e = extracted_lat_range[1]
###   check loop
        sub_domain_file_list = []
        sub_domain_lon_list = []
        sub_domain_lat_list = []
        sub_domain_utc_list = []

        for file_name in full_path_file_list:#(6,7):
            re = nc.Dataset(file_name, 'r',  format='NETCDF4_CLASSIC')
            if int(re.bad)<0.5:
                re = nc.Dataset(file_name, 'r',  format='NETCDF4_CLASSIC')
                ro_lon = re.variables['lon'][0]
                ro_lon[ro_lon<0]+=360
                ro_lat = re.variables['lat'][0]

                if ro_lon>=lon_s and ro_lon<=lon_e and ro_lat>=lat_s and ro_lat<=lat_e:
                   sub_domain_file_list.append(file_name)
                   sub_domain_lon_list.append(float(ro_lon))
                   sub_domain_lat_list.append(float(ro_lat))
        if lonlat_list:
            return(sub_domain_file_list,sub_domain_lon_list,sub_domain_lat_list)
        else:
            return(sub_domain_file_list)

    def read_profile(self,sub_domain_file_list,var):
        import netCDF4 as nc
        import numpy as np
        sub_domain_file_list = self.check_list(sub_domain_file_list)
        var = self.check_list(var)
        profiles = []
        if len(var) > 1:
           print('Var list length > 1, beware whether the output value is correct')
        for file_name  in sub_domain_file_list:
            re = nc.Dataset(file_name, 'r',  format='NETCDF4_CLASSIC')
           
            for var_name in var:
                var_prof = np.array(re.variables[var_name])
                if var_name == 'lon':
                   var_prof[var_prof<0]+=360
                if var_name == 'Temp':
                   var_prof = var_prof + 273.15
                profiles.append(var_prof)
        return(profiles)


    def ro_file_time(self,full_path_file_list,full_info=False):
        full_path_file_list = self.check_list(full_path_file_list)
        year = []
        julian = []
        utc_hour = []
        utc_min = []

        for file_name in full_path_file_list:
            split_name = file_name.split('/')[-1]
            yyear = split_name.split('.')[1]
            jjulian = split_name.split('.')[2]
            hhour = split_name.split('.')[3]
            mmin = split_name.split('.')[4]
            year.append(int(yyear))
            julian.append(int(jjulian))
            utc_hour.append(int(hhour))
            utc_min.append(int(mmin))
        if full_info:
            return(year,julian,utc_hour,utc_min)
        else:
            return(utc_hour)

    def plot_ro_distribution(self,ro_file_list,
                             rgb_array=[],plotting_lon_list=[],plotting_lat_list=[],
                             cloudsat_file=[], 
                             extracted_lon=[],extracted_lat=[],
                             extracted_hdf_date=[],extracted_hdf_id=[],extracted_hdf_granule=[],
                             extracted_region_start_time=[],
                             extracted_region_end_time=[],
                             extracted_lon_range=[],extracted_lat_range=[],
                             coast_line_color='gold',
                             lonlat_c='k',lonlat_width=0.01,lonlat_order=1,
                             lonlat_step=4,font_size=24,
                             loc='left',pad=12,
                             figure_title='RO distribution',rgb_product='',rgb_time='',
                             prefix='ro_distribution',figure_name=[],
                             dpi=300, figure_path='ro_fig',
                             sounding_station=False,
                             prof_num=True,utc_color=True,
                             granule_label=False,id_label=False,
                             date_label=False,extracted_time=False,
                             save_fig=True):
        import matplotlib.pyplot as plt
        import numpy as np    
        import matplotlib.cm as cm
        import matplotlib as mpl
        from matplotlib.colors import ListedColormap
        from ael_satellite_tools.preprocess import CloudSat
        cloudsat = CloudSat(work_path=[],
                            lat_range=self.lat_range,
                            lon_range=self.lon_range)
        ro_file_list = self.check_list(ro_file_list)
        rgb_array = self.check_list(rgb_array)
        plotting_lon_list = self.check_list(plotting_lon_list)
        plotting_lat_list = self.check_list(plotting_lat_list)

### call outer function
        extracted_file_list,extracted_lon_list,extracted_lat_list \
        = self.sub_domain_check(ro_file_list,
                                extracted_lon_range,
                                extracted_lat_range,
                                lonlat_list=True)
        extracted_utc_list = self.ro_file_time(extracted_file_list)
###
### plotting unit
### lon lat info
        lon_s,lon_e,lat_s,lat_e,lon_range,lat_range, \
        lon_step,lat_step = self.plot_lonlat_info(lonlat_step,
                                                  extracted_lon_range,
                                                  extracted_lat_range)
### basemap
        ax1,m = self.plot_basemap_unit(lon_s,lon_e,lat_s,lat_e,
                                     lon_range,lat_range,
                                     lon_step,lat_step,
                                     coast_line_color=coast_line_color,
                                     lonlat_c=lonlat_c,lonlat_width=lonlat_width,
                                     lonlat_order=lonlat_order,
                                     font_size=font_size)
### sounding sataion
        if sounding_station:
###         TAI
            new_tpe = [121.52, 24.9593]        ## 46692 
            pj_islet = [122.079744, 25.627975] ## 46695
            hua = [121.613275, 23.975128]      ## 46699
            dongsha = [116.72, 20.70000]       ## 46810 
###         JPN
            ish = [124.160, 24.330]            ## 47918
            minamidaito = [131.230, 25.830]    ## 47945
            naze = [129.550, 28.380]           ## 47909
            kago = [130.550, 31.550]           ## 47827
            fukuo = [130.380, 33.580]          ## 47807
###         KOREA
            heuksando = [125.451, 34.687]      ## 47169
            stat_c = 'fuchsia'
###
            b0 = m.scatter(new_tpe[0],new_tpe[1],210,marker='*',c=stat_c,zorder=5)
            b1 = m.scatter(pj_islet[0],pj_islet[1],210,marker='*',c=stat_c,zorder=5)
            b2 = m.scatter(hua[0],hua[1],210,marker='*',c=stat_c,zorder=5)
            b3 = m.scatter(dongsha[0],dongsha[1],210,marker='*',c=stat_c,zorder=5)
###
            b4 = m.scatter(ish[0],ish[1],210,marker='*',c=stat_c,zorder=5)
            b5 = m.scatter(minamidaito[0],minamidaito[1],210,marker='*',c=stat_c,zorder=5)
            b6 = m.scatter(naze[0],naze[1],210,marker='*',c=stat_c,zorder=5)
            b7 = m.scatter(kago[0],kago[1],210,marker='*',c=stat_c,zorder=5)
            b8 = m.scatter(fukuo[0],fukuo[1],210,marker='*',c=stat_c,zorder=5)
###
            b9 = m.scatter(heuksando[0],heuksando[1],210,marker='*',c=stat_c,zorder=5)

### himawari 
        if len(rgb_array)>0:
            fig_num = 0
            xi, yi = np.meshgrid(plotting_lon_list[0],plotting_lat_list[0])
            rgb_array_shape = rgb_array[fig_num].shape
            if len(rgb_array_shape)>2:
                m.pcolormesh(xi,yi,rgb_array[fig_num])
            else: 
                bt_min = np.min(rgb_array[fig_num])+10
                bt_max = np.max(rgb_array[fig_num])-10
                m.pcolormesh(xi,yi,rgb_array[fig_num],vmin=bt_min,vmax=bt_max,cmap='Greys')

### cloudsat
        if len(cloudsat_file)>0:
            cloudsat_file_list,region_mask_list \
            = cloudsat.sub_domain_check(cloudsat_file,mask_output=True)

            plot_lon,plot_lat, hdf_date,hdf_id,hdf_granule,\
            region_start_time,region_end_time \
            = cloudsat.read_geometric_info(cloudsat_file_list,
                                           region_mask_list,
                                           lonlat_info=True)
### plot track main
            cloudsat.plot_track_unit(ax1,plot_lon,plot_lat,
                                     hdf_date,hdf_id,hdf_granule,
                                     lon_range,lat_range,
                                     font_size=font_size,granule_label=granule_label,
                                     id_label=id_label,date_label=date_label)
        if len(extracted_lon) > 0:
            cloudsat.plot_track_unit(ax1,extracted_lon,extracted_lat,
                                 extracted_hdf_date,extracted_hdf_id,extracted_hdf_granule,
                                 lon_range,lat_range,line_color='mediumblue',
                                 font_size=font_size,granule_label=False,
                                 id_label=False,date_label=False)


###
        if utc_color:
            cmap1 = mpl.cm.jet(np.linspace(0,1,24))
            #cmap2 = mpl.cm.nipy_spectral(np.linspace(0,1,24))
            blue1 = [0.22,0.39,0.90,1.]
            green = [0.5,1.,0.,1.]
            red = [0.92, 0., 0., 1.]
            med_orchid = [0.73,0.33,0.83,1.]
            self_cmap = np.vstack((blue1,cmap1[8],green,cmap1[17],red,med_orchid))
            cmap = mpl.colors.ListedColormap(self_cmap)
            cmap.set_extremes(over=cmap1[23])
            ro_plot = m.scatter(extracted_lon_list,extracted_lat_list,
                                140,c=extracted_utc_list,cmap=cmap,
                                vmin=0,vmax=24,zorder=4)
        else:
            ro_plot = m.scatter(extracted_lon_list,extracted_lat_list,
                                140,c='r',
                                zorder=4)

        if prof_num:
## profile number plotting
            extracted_lon_list = np.array(extracted_lon_list)
            extracted_lat_list = np.array(extracted_lat_list)
 
            for i in range(0,len(extracted_file_list)):#len(ffiles)):
                clos_lon = -0.8
                clos_lat = 0.8
                scale = 1/(0.3*(np.sqrt(lon_range*lon_range + lat_range*lat_range)/18))
                if len(extracted_file_list) > 1:
                    dis_lon = np.abs(extracted_lon_list-extracted_lon_list[i])
                    dis_lat = np.abs(extracted_lat_list-extracted_lat_list[i])
                    total_dis = np.sqrt(dis_lon*dis_lon + dis_lat*dis_lat)
                    clos_dis = np.sort(total_dis)[1]
                    scale = clos_dis/(0.3*(np.sqrt(lon_range*lon_range + lat_range*lat_range)/18))
                    clos_indx = np.argsort(total_dis)[1]
                    clos_lon = extracted_lon_list[clos_indx] - extracted_lon_list[i]
                    clos_lat = extracted_lat_list[clos_indx] - extracted_lat_list[i]
                #print(clos_lon,clos_lat)
                #print(extracted_lon_list[i],extracted_lat_list[i])
                
                if len(rgb_array)>0:
                    text_color = 'lime'
                else:
                    text_color = 'seagreen'
                plt.text(extracted_lon_list[i]-(clos_lon/scale),
                         extracted_lat_list[i]-(clos_lat/scale),
                         str(i),c=text_color,fontsize=20,fontweight='bold',zorder=4,
                         horizontalalignment='center', verticalalignment='center')
###
### plot title
        ro_title = ''
        cloudsat_title = '' 
        hima_title = ''
        if len(self.date_range) > 0:
            ro_title = ''+figure_title+'' +' '+self.date_range[0][0:5] + self.date_range[0][5:] +'~'+ self.date_range[1][5:]
        else:
            ro_title = ''+figure_title+''
        if len(cloudsat_file) > 0:
            cloudsat_title = 'CloudSat track '+ \
                             str(region_start_time[0][0]) + ':' + \
                             str(region_start_time[0][1]) + 'UTC' + '~' + \
                             str(region_end_time[0][0]) + ':' + \
                             str(region_end_time[0][1]) + 'UTC'
            
        if len(rgb_array) > 0:
            hima_title = rgb_product+' '+rgb_time
#        title_date =  self.date_range[0] +'~'+ self.date_range[1]
#        title_date = ''
        title = ''+ro_title+'' 
        if len(cloudsat_file) > 0:
             title = ''+title+' \n'+cloudsat_title+''
        if len(rgb_array) > 0: 
             title = ''+title+' \n'+hima_title+''
        self.plot_title_unit(title,loc=loc,pad=pad,font_size=font_size-2)
### colorBar
        if utc_color:
            if lon_range > lat_range:
                hei_size = 13.3/13.3
            else:
                hei_size = 13.4*(lat_range/lon_range)

        #cbar_ax = ax1.add_axes([0.82, 1-0.18*(hei_size), 0.15, 0.03*(hei_size)])
            c_ref = plt.colorbar(ro_plot,ticks=[0, 4, 8, 12, 16, 20, 24],
                                 pad=0.01,fraction=0.012,aspect=20,anchor=(0.95,0.06),
                                 orientation='horizontal', location='top')
            c_ref.ax.tick_params(labelsize=14)
            c_ref.set_label(label='UTC hour', size=12, weight='bold')

### save fig
        if len(self.date_range) > 0:        
            save_date1 = self.date_range[0][0:-3].replace('/','').replace(' ','').lower()
            save_date2 = self.date_range[1][0:-3].replace('/','').replace(' ','').lower()
            figure_main = save_date1+'_'+save_date2
        else:
            figure_main = ''
        self.save_fig_unit(figure_main,figure_name=figure_name,
                           prefix=prefix,
                           dpi=dpi, figure_path=figure_path,save_fig=save_fig)

    def plot_ro_profile(self,ro_file_list,
                        prof_num = [],
                        moist_type='sph',height_type='MSL_alt',
                        height_range=[],height_step=5,
                        fontsize=16,
                        prefix='ro_profile',figure_name=[],
                        dpi=300, figure_path='ro_fig',
                        profile_label=False,save_fig=True):
        import matplotlib.pyplot as plt
        import numpy as np
        #temp_prof = self.check_list(temp_prof)
        #moist_prof = self.check_list(moist_prof)
        #vertical_prof = self.check_list(vertical_prof)
        ro_file_list = self.check_list(ro_file_list)
        prof_num = self.check_list(prof_num)
               
        moist1,moist2,moist_xlabel,moist_title = self.profile_setting(moist_type)
        temp1,temp2,temp_xlabel,temp_title = self.profile_setting('Temp')
### decide target profiles
        if len(prof_num)<1:
            read_s = 0
            read_e = len(ro_file_list)
            read_list = np.arange(read_s,read_e)
        if len(prof_num)>=1:
            read_list = prof_num 
### read data
        for i in read_list:
            temp_prof = self.read_profile(ro_file_list[i],'Temp')
            moist_prof = self.read_profile(ro_file_list[i],moist_type)
            vertical_prof = self.read_profile(ro_file_list[i],height_type)
            lat = self.read_profile(ro_file_list[i],'lat')[0][0]
            lon = self.read_profile(ro_file_list[i],'lon')[0][0]
 
### some plotting info prepare
            year,julian,utc_hour,utc_min = self.ro_file_time(ro_file_list[i],full_info=True)
            date = self.julian_to_monthdate(year,julian)
#
            ver1,ver2,ylabel,hei_ticks = \
            self.vertical_axis_setting(vertical_prof,height_range,height_step)

### plot figure 
            fig,ax1 = plt.subplots(1,2,figsize=(10,6.8),tight_layout=True)
### subplot Temp.
            ax1[0].plot(temp_prof[0],vertical_prof[0],linewidth=4,label='RO')
            ax1[0].axis([temp1,temp2,ver1,ver2])
            if ver1 - ver2 < 0:
                ax1[0].set_yticks(hei_ticks)
            ax1[0].grid(axis='y')
            ax1[0].tick_params(axis='both', labelsize=fontsize)
            ax1[0].set_ylabel(ylabel,fontsize=fontsize)
            ax1[0].set_xlabel(temp_xlabel,fontsize=fontsize)
            if profile_label:
                ax1[0].legend(fontsize=fontsize-2)
### 
            str_hour = str(utc_hour[0])
            str_min = str(utc_min[0])
            str_prof_num = str(i)
            if utc_hour[0] < 10:
                str_hour = '0'+str_hour
            if utc_min[0] < 10:
                str_min = '0'+str_min
            if i < 10 :
                str_prof_num = '0'+str_prof_num
            title_time = date[0] +' '+ str_hour+':'+str_min + 'UTC'
            pos = u''+str(np.around(lat,1))+'\u00B0N, ' +str(np.around(lon,1)) +'\u00B0E'
###
            ax1[0].set_title('RO Prof.  (' +str(i)+ ')  @'+pos+' \n'+temp_title+'',loc='left',fontsize=fontsize+2)
### subplot moisture profile
            ax1[1].plot(moist_prof[0],vertical_prof[0],linewidth=4,label='RO')
            ax1[1].axis([moist1,moist2,ver1,ver2])
            if ver1 - ver2 < 0:
                ax1[1].set_yticks(hei_ticks)
            ax1[1].grid(axis='y')
            ax1[1].tick_params(axis='both', labelsize=fontsize)
            ax1[1].set_ylabel(ylabel,fontsize=fontsize)
            ax1[1].set_xlabel(moist_xlabel,fontsize=fontsize)
            if profile_label:
                ax1[1].legend(fontsize=fontsize-2)
            ax1[1].set_title(''+title_time+'\n'+moist_title+' ',loc='left',fontsize=fontsize+2) 
### save figure
            figure_main = str_prof_num+'_'+date[0].replace('/','')+'_'+str_hour+str_min
            self.save_fig_unit(figure_main,figure_name=figure_name,
                               prefix=prefix,
                               dpi=dpi, figure_path=figure_path,save_fig=save_fig)
    def plot_profile_unit(self,ax1,profile,vertical_prof,var_type='Temp',label='RO',
                          height_range=[0,20],height_step=6,
                          sub_plot_num=2,plot_num=1,
                          linewidth=4,alpha=1,fontsize=16,
                          profile_label=False,axis_setting=True):
        var1,var2,var_xlabel,var_title = self.profile_setting(var_type)
        profile = self.check_list(profile)
        vertical_prof = self.check_list(vertical_prof)
        if sub_plot_num < 2:
            ax1 = [ax1]
### subplot Temp.
        for i in range(0,len(profile)):
            ax1[plot_num-1].plot(profile[i],vertical_prof[i],linewidth=linewidth,label=label,alpha=alpha)
        if axis_setting:
            ver1,ver2,ylabel,hei_ticks = \
                self.vertical_axis_setting(vertical_prof[i],
                                           height_range=height_range,
                                           height_step=height_step)
            ax1[plot_num-1].axis([var1,var2,ver1,ver2])
            if ver1 - ver2 < 0:
                ax1[plot_num-1].set_yticks(hei_ticks)
            ax1[plot_num-1].set_ylabel(ylabel,fontsize=fontsize)
            ax1[plot_num-1].grid(axis='y')
            ax1[plot_num-1].tick_params(axis='both', labelsize=fontsize)
            ax1[plot_num-1].set_xlabel(var_xlabel,fontsize=fontsize)
        if profile_label:
            ax1[plot_num-1].legend(fontsize=fontsize-2)

    def vertical_axis_setting(self,vertical_prof,height_range=[],height_step=5):
        import numpy as np
        vertical_prof = self.check_list(vertical_prof)
        hei1 = vertical_prof[0][0]
        hei2 = vertical_prof[0][-1]
        if len(height_range) < 1:
            if hei1 - hei2 > 0:
                ver1 = 1030
                ver2 = 5
            else:
                ver1 = 0
                ver2 = 20
        else:
            ver1 = height_range[0]
            ver2 = height_range[1]
        if ver1 - ver2 > 0:
            ylabel = '[hPa]'
            hei_jump = int(np.floor((ver2 - ver1)/(height_step-1)))
            hei_range = list(np.arange(ver1,ver2+1,hei_jump))
            hei_ticks = hei_range[:]
        else:
            ylabel = '[km]'
            hei_jump = int(np.floor((ver2 - ver1)/(height_step-1)))
            if hei_jump == 0:
                hei_jump = 1
            hei_range = list(np.arange(ver1,ver2+1,hei_jump))
            hei_ticks = hei_range[:]
            if len(height_range) > 0:
                if height_range[0] < 1:
                     hei_ticks = [1]
                     hei_ticks.extend(hei_range[1:])
        return(ver1,ver2,ylabel,hei_ticks)

    def profile_setting(self,var_type):
        if var_type == 'Temp':
            moist1 = 185
            moist2 = 310
            moist_xlabel = 'Temperature [K]'
            moist_title = 'Temperature'
        if var_type == 'sph':
            moist1 = -1
            moist2 = 20
            moist_xlabel = 'q [g/kg]'
            moist_title = 'q'
        if var_type == 'Vp':
            moist1 = -1
            moist2 = 30
            moist_xlabel = 'Vapor Pressure [hPa]'
            moist_title = 'VP'
        if var_type == 'rh':
            moist1 = -2
            moist2 = 102
            moist_xlabel = 'RH [%]'
            moist_title = 'RH'
        return(moist1,moist2,moist_xlabel,moist_title)

    def julian_to_monthdate(self,year,julian):
        from datetime import datetime, timedelta
        year = self.check_list(year)
        julian = self.check_list(julian)
        date = []
        for i in range(0,len(year)):
            start_date = datetime(year[i], 1, 1)
            date_result = start_date + timedelta(days=julian[i] - 1)
            date.append(date_result.strftime("%Y/%m/%d"))
        return(date)

    def ro_information(self,detail=False):
        print('Time: UTC time')
        print('Variable unit:')
        print('Geometric var')
        print('----------------------------')
        print('lat: -90 ~ 90 ')
        print('lon: -180 ~ 180 (module output adjust to 0 ~ 360)')
        print('----------------------------')
        print('Altitude var')
        print('----------------------------')
        print('MSL_alt (height above MSL): km')
        print('gph (geopotential height): km')
        print('Pres (pressure): hPa')
        print('----------------------------')
        print('Thermo var')
        print('----------------------------')
        print(u'Temp (temperature): \u00B0C (module output adjust to K)')
        print('sph (specific humidity): g/kg')
        print('Vp (vapor pressure): mb')
        print('rh (relative humidity):  %')
        print('----------------------------')
        print('missing value: -999.f')
        if detail:
            print('global attributes: bad = "0" profile generated')
            print('global attributes: bad = "1" no profile generated')
            print('----------------------------')
            print('Extra var')
            print('----------------------------')
            print('ref (Refractivity):  N')
            print(u'temp_dry (DRY T derived from optimized BA):  \u00B0C')
            print('pres_dry (DRY p derived from optimized BA):  hPa')

