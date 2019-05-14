function PRE_proc_2bats

clear all
close all

eval ('params_in')


day_rows = 1:18;
PRE_create_day_structs (p_in,day_rows)

cell_rows = 1:119;
PRE_create_cells_struct (p_in,cell_rows)

end