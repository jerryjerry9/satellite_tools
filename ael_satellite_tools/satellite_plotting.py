"""
satellite_plotting.py

This module provides plotting functions for satellite associated modules

"""

class Plotting_tools:
    def __init__(self,satellite=[]):
        self.satellite = satellite

    def plot_lonlat_info(self,lonlat_step,extracted_lon_range=[],extracted_lat_range=[]):
        import numpy as np
        if len(extracted_lon_range)<1:
            lon_s = self.lon_range[0]
            lon_e = self.lon_range[1]
        else:
            lon_s = extracted_lon_range[0]
            lon_e = extracted_lon_range[1]

        if  len(extracted_lat_range)<1:
            lat_s = self.lat_range[0]
            lat_e = self.lat_range[1]
        else:
            lat_s = extracted_lat_range[0]
            lat_e = extracted_lat_range[1]

        lon_range = lon_e - lon_s
        lat_range = lat_e - lat_s
        lon_step = np.floor((lon_range)/(lonlat_step-1))
        lat_step = np.floor((lat_range)/(lonlat_step-1))
        return(lon_s,lon_e,lat_s,lat_e,lon_range,lat_range,lon_step,lat_step)

    def plot_basemap_unit(self,lon_s,lon_e,lat_s,lat_e,lon_range,lat_range,lon_step,lat_step,coast_line_color='olive',lonlat_c='k',lonlat_width=0.01,lonlat_order=1,font_size=24):
        import numpy as np
        import matplotlib.pyplot as plt
        from mpl_toolkits.basemap import Basemap
        if lon_range > lat_range:
            fig,ax1 = plt.subplots(1,1,figsize=(13*(lon_range/lat_range),13.3),tight_layout=True)
        else:
            fig,ax1 = plt.subplots(1,1,figsize=(13,13.4*(lat_range/lon_range)),tight_layout=True)
        m = Basemap(llcrnrlon=lon_s, urcrnrlon=lon_e, llcrnrlat=lat_s, urcrnrlat=lat_e,resolution='i')
        m.drawparallels(np.arange(lat_s, lat_e+1, lat_step),
                        labels=[1, 0, 0, 0], linewidth=lonlat_width, color=lonlat_c,
                        fontsize=font_size-2,zorder=lonlat_order)
        m.drawmeridians(np.arange(lon_s, lon_e+1, lon_step),
                        labels=[0, 0, 0, 1], linewidth=lonlat_width, color=lonlat_c,
                        fontsize=font_size-2,zorder=lonlat_order)
        m.drawcoastlines(linewidth=1.3,color=coast_line_color,zorder=3)
        return(ax1,m)
   
    def plot_profile_unit(self,profile_data,x_axis,reference_hei,shading='contourf',color_map=[],value_range=[0,100],extend='min',zorder=2):
        import numpy as np
        import matplotlib.pyplot as plt
        if shading == 'pcolormesh':
            prof = plt.pcolormesh(x_axis,reference_hei,profile_data,
                                cmap=color_map,
                                vmin=value_range[0],vmax=value_range[1],
                                zorder=2)
        if shading == 'contourf':
            prof = plt.contourf(x_axis,reference_hei,profile_data,
                                cmap=color_map,extend=extend,
                                levels=np.arange(value_range[0],value_range[1]+0.05,0.1),zorder=2)
        return(prof)

    def plot_title_unit(self,title,loc='left',pad=40,font_size=24):
### figure info
        import matplotlib.pyplot as plt
        plt.title(title,loc=loc,fontsize=font_size,pad=pad)

    def save_fig_unit(self,figure_main,figure_name=[],prefix=[],
                      dpi=300, figure_path=[],save_fig=True):
        from pathlib import Path
        import matplotlib.pyplot as plt
        if save_fig:
            if figure_path == None or figure_path == []:
                figure_path = ''+str(self.work_path)+'/'
            else:
                Path(figure_path).mkdir(parents=True, exist_ok=True)
                figure_path = ''+str(figure_path)+'/'

            if len(figure_name)>0:
                figure_main = ''+figure_name+'_'+figure_main+''
            plt.savefig(''+figure_path+''+prefix+'_'+figure_main+'.png',dpi=dpi)

    def plot_timezone_unit(self,ax1,lon_s,lon_e,lat_s,lat_e,lat_range,font_size=24,utc_label=True):
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.patches import Polygon
        import matplotlib.cm as cm
### time zone
        lon_min = np.arange(-7.5,338,15)
        lon_max = np.arange(7.5,353,15)
        lat_min, lat_max = lat_s, lat_e
        color_step, local_ref_time = self.cloudsat_string_info()
        cmap = cm.get_cmap('Blues_r', len(lon_min)*3)
        zone_num = int(lon_e/15)+1
        if zone_num >24:
           zone_num = 24
        for j in range(int(lon_s/15),zone_num):
            polygon = [(lon_min[j], -90), (lon_max[j], -90), (lon_max[j], 90), (lon_min[j], 90)]
            color = cmap(color_step[j]*3 / len(lon_min))
            poly_map = Polygon(polygon, edgecolor='thistle', facecolor=color, alpha=1, linewidth=2)
            ax1.add_patch(poly_map)
            utc_pos = (lon_max[j] + lon_min[j])/2
            if utc_pos < lon_s:
                utc_pos = lon_s
            if utc_pos > lon_e:
                utc_pos = lon_e
            if utc_label:
                plt.text(utc_pos,lat_max+(lat_range/180),''+local_ref_time[j]+'UTC',c='teal',
                         fontsize=font_size-2, fontweight='bold',
                         va='bottom',ha='center',
                         bbox=dict(facecolor='ghostwhite', edgecolor='indigo', boxstyle='round,pad=0.1'),zorder=6)

    """
    color map
    """
    def customized_colormap(self,color_map_name):
        color_map_name = self.check_list(color_map_name)
        color_map = color_map_name[0] 
        method = getattr(self, color_map)
        output_color_map = method()
        return(output_color_map)

    def light_YlGnBu(self):
        import matplotlib.colors as mcolors
        ### colormap setting
        light_YlGnBu = mcolors.LinearSegmentedColormap.from_list(
            "light_YlGnBu", ["#ffffcc", "#a1dab4", "#41b6c4", "#2c7fb8", "#2870c9"], N=256)
        light_YlGnBu.set_under((1, 1, 1, 0))
        return(light_YlGnBu)


       
