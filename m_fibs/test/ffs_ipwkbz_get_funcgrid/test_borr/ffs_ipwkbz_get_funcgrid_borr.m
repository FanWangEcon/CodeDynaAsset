%% Generate Choice Grids for Borrowing for IPKWZ

close all

%% Choice Grids, Zoom in to Borrow 
% this is the default grid

% get param_map and support_map
it_param_set = 4;
[param_map, support_map] = ffs_fibs_set_default_param(it_param_set);
support_map('bl_graph_funcgrids') = true;
support_map('bl_display_funcgrids') = false;

% to be able to visually see choice grid points
param_map('fl_b_bd') = -20;
param_map('fl_w_min') = param_map('fl_b_bd');
param_map('it_w_perc_n') = 5;
param_map('it_ak_perc_n') = 5;
param_map('fl_w_interp_grid_gap') = 3;
param_map('fl_coh_interp_grid_gap') = 3;

% max
param_map('fl_w_max') = 0;
param_map('fl_k_max') = (param_map('fl_w_max') - param_map('fl_b_bd'));

% minimum income
param_map('fl_w') = 0.4;

% set very large borrowing interest rate, differential borr save rates
param_map('fl_r_save') = 0.025;
param_map('fl_r_borr') = 0.035;

% call program
ffs_fibs_get_funcgrid(param_map, support_map)

%% Choice Grids, Low Borrow Interest, 35 percent
% this is the default grid

% get param_map and support_map
it_param_set = 4;
[param_map, support_map] = ffs_fibs_set_default_param(it_param_set);
support_map('bl_graph_funcgrids') = true;
support_map('bl_display_funcgrids') = false;

% to be able to visually see choice grid points
param_map('fl_b_bd') = -20;
param_map('fl_w_min') = param_map('fl_b_bd');
param_map('it_w_perc_n') = 5;
param_map('it_ak_perc_n') = 5;
param_map('fl_w_interp_grid_gap') = 3;
param_map('fl_coh_interp_grid_gap') = 3;

% minimum income
param_map('fl_w') = 0.4;

% set very large borrowing interest rate, differential borr save rates
param_map('fl_r_save') = 0.025;
param_map('fl_r_borr') = 0.035;

% call program
ffs_fibs_get_funcgrid(param_map, support_map)

%% Choice Grids, High Borrow Interest, 35 percent
% this is the default grid

% get param_map and support_map
it_param_set = 4;
[param_map, support_map] = ffs_fibs_set_default_param(it_param_set);
support_map('bl_graph_funcgrids') = true;
support_map('bl_display_funcgrids') = false;

% to be able to visually see choice grid points
param_map('fl_b_bd') = -20;
param_map('fl_w_min') = param_map('fl_b_bd');
param_map('it_w_perc_n') = 5;
param_map('it_ak_perc_n') = 5;
param_map('fl_w_interp_grid_gap') = 3;
param_map('fl_coh_interp_grid_gap') = 3;

% minimum income
param_map('fl_w') = 0.4;

% set very large borrowing interest rate, differential borr save rates
param_map('fl_r_save') = 0.025;
param_map('fl_r_borr') = 0.35;

% call program
ffs_fibs_get_funcgrid(param_map, support_map)

%% Choice Grids, High Borrow Interest, 35 percent, Low Borrow Bound
% this is the default grid

% get param_map and support_map
it_param_set = 4;
[param_map, support_map] = ffs_fibs_set_default_param(it_param_set);
support_map('bl_graph_funcgrids') = true;
support_map('bl_display_funcgrids') = false;

% to be able to visually see choice grid points
param_map('fl_b_bd') = -100;
param_map('fl_w_min') = param_map('fl_b_bd');
param_map('it_w_perc_n') = 5;
param_map('it_ak_perc_n') = 5;
param_map('fl_w_interp_grid_gap') = 3;
param_map('fl_coh_interp_grid_gap') = 3;

% minimum income
param_map('fl_w') = 0.4;

% set very large borrowing interest rate, differential borr save rates
param_map('fl_r_save') = 0.025;
param_map('fl_r_borr') = 0.35;

% call program
ffs_fibs_get_funcgrid(param_map, support_map)