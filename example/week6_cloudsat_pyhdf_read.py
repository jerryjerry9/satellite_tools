from pyhdf.HDF import HDF, HC
from pyhdf import V
from pyhdf.SD  import SD, SDC

geoprof_file = ['/data/dadm1/obs/CloudSat/hdf-GEOPROF/2016/2016359034029_56694_CS_2B-GEOPROF_GRANULE_P1_R05_E06_F01.hdf']

### read Vector data variable list
reading_vdata = HDF(geoprof_file[0], HC.READ).vstart()
vdata_list = reading_vdata.vdatainfo()

print(vdata_list)

### read Vector data 
reading_vdata = HDF(geoprof_file[0], HC.READ).vstart()
vdata_list = reading_vdata.vdatainfo()
target_vdata = 'Latitude'
for ref in vdata_list:
    vdata_name = ref[0]
    if vdata_name == target_vdata:
        vdata_length = ref[3]
        Var_data = reading_vdata.attach(target_vdata)
        var_data = Var_data.read(vdata_length)
        Var_data.detach()
        break
reading_vdata.end()

### read Scientific Data variable list
hdfFile = SD(geoprof_file[0], SDC.READ)
dsets = hdfFile.datasets()

print(dsets)

### read Scientific Data
hdfFile = SD(geoprof_file[0], SDC.READ)
Var_data = hdfFile.select('Radar_Reflectivity')
var_data = Var_data[:]
hdfFile.end
