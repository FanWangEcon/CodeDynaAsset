%% Test Borrowing *Default* (Risky + Safe Asset + Save + Borr + R Shock + Interpolated-Percentage), Cross Test
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%

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

%% Medium CROSS Simulation Inf Preference (Limited Graphs)
% Generate Arrays of Parameter Values to Loop Over

cl_st_param_keys = {'fl_crra', 'fl_beta'};

it_simu_vec_len = 15;
param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map('fl_crra') = linspace(1, 5, it_simu_vec_len);
param_tstar_map('fl_beta') = linspace(0.87, 0.97, it_simu_vec_len);

it_size_type = 2;
ar_it_plot_sets = [1001, 1002, 1003];
bl_simu_cross = 'c';

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);

%% Medium CROSS Simulation For R Borr Save (Limited Graphs)
% Generate Arrays of Parameter Values to Loop Over

cl_st_param_keys = {'fl_r_fsv', 'fl_r_fbr'};

it_simu_vec_len = 15;
param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map('fl_r_fsv') = linspace(0.00,    0.05, it_simu_vec_len);
param_tstar_map('fl_r_fbr') = linspace(0.03, 0.08, it_simu_vec_len);

it_size_type = 2;
ar_it_plot_sets = [1001, 1002, 1003];
bl_simu_cross = 'c';

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);

%% Medium CROSS Simulation Inf Average Rate and Borrow Gap (Limited Graphs)

cl_st_param_keys = {'fl_z_r_infbr_poiss_mean', 'fl_forbrblk_gap'};

it_simu_vec_len = 15;
param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map('fl_z_r_infbr_poiss_mean') = linspace(2, 10, it_simu_vec_len);
param_tstar_map('fl_forbrblk_gap') = linspace(-0.1, -2.5, it_simu_vec_len);

it_size_type = 2;
ar_it_plot_sets = [1001, 1002, 1003];
bl_simu_cross = 'c';

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);
