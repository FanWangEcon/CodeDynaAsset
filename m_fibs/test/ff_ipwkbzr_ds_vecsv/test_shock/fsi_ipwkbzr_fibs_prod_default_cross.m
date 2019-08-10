%% Test Preference *Default* (Risky + Safe Asset + Save + Borr + R Shock + Interpolated-Percentage), Test over Discount and Risk-Aversion Arrays (cross)

%% Set Shared Parameters
close all;
clear all;

% Borrowing/Savings Parameters
bl_default = true;

%% Simulate and Graph
% Note: as for example _fl_beta_ increases, willingness to save increases,
% leading to higher savings, which will exceed the benchmark max grid
% point. So to allow for higher beta, dramatically higher max savings bound
% is needed.

% Set which to graph, simulate over which variables
cl_st_param_keys = {'fl_alpha', 'fl_delta'};

% Generate Benchmark Parameters
it_param_set = 9;
[param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);

% Borrowing Parameters
param_map('bl_default') = bl_default;

% Support Parameters
support_map('bl_replacefile') = false;
support_map('bl_graph_onebyones') = true;
support_map('bl_display_graph_stats') = false;
support_map('st_mat_test_prefix') = ['dft_'];

% Generate Arrays of Parameter Values to Loop Over
it_simu_vec_len = 15;
param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map('fl_alpha') = linspace(0.30, 0.50, it_simu_vec_len);
param_tstar_map('fl_delta') = linspace(0.02, 0.14, it_simu_vec_len);

%% Quick CROSS Simulation (Limited Graphs)
% it_size_type = 1;
% ar_it_plot_sets = [3,4,102, 152,104,106];
% bl_simu_cross = 'c';
% 
% % Simulate along parameters
% ff_az_test_analyze( ...
%     ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
%     param_map, support_map, param_tstar_map);
% 
% close all;

%% Medium CROSS Simulation (Limited Graphs)
it_size_type = 2;
ar_it_plot_sets = [3,4,102, 152,104,106, 10, 201, 1001, 1002, 1003];
bl_simu_cross = 'c';

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);

close all;

%% Medium GRID Simulation (GRID Limited Graphs)
% it_size_type = 3;
% ar_it_plot_sets = [3,4,102, 152,104,106];
% bl_simu_cross = 'g';
% 
% % Simulate along parameters
% ff_az_test_analyze( ...
%     ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
%     param_map, support_map, param_tstar_map);
% 
% close all

%% Denser CROSS Simulation
% it_size_type = 3;
% ar_it_plot_sets = [51,52,53,54, 5,6,103,153, 61,62,63,64];
% bl_simu_cross = 'c';
% 
% % Simulate along parameters
% ff_az_test_analyze( ...
%     ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
%     param_map, support_map, param_tstar_map);
% 
% close all;

%% Denser GRID Simulation
% it_size_type = 3;
% ar_it_plot_sets = [1,2,101,151, 3,4,102,152, 5,6,103,153, 51,52,53,54, 201,205,207,209, 104,105,106,10];
% bl_simu_cross = 'g';
% 
% % Simulate along parameters
% ff_az_test_analyze( ...
%     ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
%     param_map, support_map, param_tstar_map);
% 
% close all