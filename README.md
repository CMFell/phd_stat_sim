# phd_stat_sim

This contains two sets of code, one for simulations to test multinomial GAMs (simulation_code) and one using multinomial GAMS for modelling the Greater Fish River Canyon survey area (gfrc_modelling). 

In simulation code the files are named part1 to part 13 and should be run in order, at the top of each file you need to set a local base directory to save output and read input from previous files. This is base_dir_on_comp and is usally the first non comment line in the file.

In gfrc_modelling any of the files can be run individually. All files need a local base directory set to save output, this is base_dir_on_comp and should be the first non comment line in the file. They also assume the data for the GFRC survey is saved relative to the base directory at data_files/final_image_data.csv . this file can be downloaded from *add link*
