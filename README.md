# ael_satellite_tools

A Python module for downloading and analysis the satellite.

0.0.13
Himawari Data Storage Directory Update

Introduced a new data_end_path configuration option to allow users to control the data storage structure.

Default setting: data_end_path='day'
Data will now be stored in the format:
"data_path/compressed_data/[YYYY]/[MM]/[DD]/*.bz2"
"data_path/sub_domain_data/[YYYY]/[MM]/[DD]/*.nc"


Legacy compatibility:
Setting data_end_path='month' will retain the previous structure:
"data_path/compressed_data/[YYYY]/[MM]/*.bz2"
"data_path/sub_domain_data/[YYYY]/[MM]/*.nc"


##Usage examples:
    from ael_satellite_tools.preprocess import Himawari as Himawari
    from ael_satellite_tools.plotting import Himawari as Hima_plot

    'Current default behavior (daily folder structure)'
    himawari = Himawari(data_end_path='day')
    hima_plot = Hima_plot(data_end_path='day')

    'Previous behavior (monthly folder structure)'
    himawari = Himawari(data_end_path='month')
    hima_plot = Hima_plot(data_end_path='month')


Feature Added: himawari.move_data(current_folder='month', target_folder='day')

Introduced a helper function to move downloaded data from the month folder to the day folder.

To move data in the opposite direction (from day to month), simply reverse the keyword arguments.


Upcoming update 0.0.14

EarthCARE process module


## Installation

```bash
pip install ael_satellite_tools

