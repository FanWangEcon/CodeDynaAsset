%% Generate Choice Grids for Borrowing for IPKWZ

close all

% %% Generate Grid for Savings
% % this is the default grid
% 
% % get param_map and support_map
% it_param_set = 4;
% [param_map, support_map] = ffs_ipwkz_set_default_param(it_param_set);
% support_map('bl_graph_funcgrids') = true;
% support_map('bl_display_funcgrids') = false;
% 
% % to be able to visually see choice grid points
% param_map('fl_b_bd') = 0;
% param_map('fl_w_min') = param_map('fl_b_bd');
% param_map('it_w_perc_n') = 25;
% param_map('it_ak_perc_n') = 45;
% param_map('fl_w_interp_grid_gap') = 2;
% param_map('fl_coh_interp_grid_gap') = 2;
% 
% ffs_ipwkz_get_funcgrid(param_map, support_map);

%% Conditional on Lowest Shock, COH Maximizing K
% what is the w=k'+b' level given which, for all shocks, there exists at
% least one k grid points where the the resulting COH is above borrowing
% bound. 

%% Generate Grid for Savings
% this is the default grid

% get param_map and support_map
it_param_set = 4;
[param_map, support_map] = ffs_ipwkbz_set_default_param(it_param_set);
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
param_map('fl_r_borr') = 0.85;

% call program
ffs_ipwkbz_get_funcgrid(param_map, support_map)