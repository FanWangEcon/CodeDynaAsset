%% Test Borrowing *With Default* (Save + Borr Distribution), Grid Test
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% @seealso
%
% * test speed: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_speed/html/fsi_abz_ds_vecsv_speed.html fsi_abz_ds_vecsv_speed>
% * test joint *RANDOM*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_joint/html/fsi_abz_ds_vecsv_joint_rand.html fsi_abz_ds_vecsv_joint_rand>
%
% * test price no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_price/html/fsi_abz_ds_vecsv_price_nbc_cross.html fsi_abz_ds_vecsv_price_nbc_cross>
% * test price default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_price/html/fsi_abz_ds_vecsv_price_default_cross.html fsi_abz_ds_vecsv_price_default_cross>
%
% * test interest rate no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc.html fsi_abz_ds_vecsv_nbc>
% * test interest rate no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc_cross.html fsi_abz_ds_vecsv_nbc_cross>
% * test interest rate no default *GRID*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc_grid.html fsi_abz_ds_vecsv_nbc_grid>
% * test interest rate default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default.html fsi_abz_ds_vecsv_default>
% * test interest rate default *V(a,z)* Comparison: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default_compaz.html fsi_abz_ds_vecsv_default_compaz>
% * test interest rate default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default_cross.html fsi_abz_ds_vecsv_default_cross>
% * test interest rate default *GRID*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default_grid.html fsi_abz_ds_vecsv_default_grid>
%
% * test shock default (very low cmin): <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_default_lowcmin.html fsi_abz_ds_vecsv_shk_default_lowcmin>
% * test shock no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_nbc.html fsi_abz_ds_vecsv_shk_nbc>
% * test shock no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_nbc_cross.html fsi_abz_ds_vecsv_shk_nbc_cross>
% * test shock default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_default.html fsi_abz_ds_vecsv_shk_default>
% * test shock default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_default_cross.html fsi_abz_ds_vecsv_shk_default_cross>
%
% * test preference no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_nbc.html fsi_abz_ds_vecsv_pref_nbc>
% * test preference no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_nbc_cross.html fsi_abz_ds_vecsv_pref_nbc_cross>
% * test preference default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_default.html fsi_abz_ds_vecsv_pref_default>
% * test preference default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_default_cross.html fsi_abz_ds_vecsv_pref_default_cross>
% * test preference default (very low cmin): <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_default_lowcmin.html fsi_abz_ds_vecsv_pref_default_lowcmin>
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
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% Borrowing Parameters
param_map('bl_default') = bl_default;

% Support Parameters
support_map('bl_replacefile') = false;
support_map('bl_graph_onebyones') = true;
support_map('bl_display_graph_stats') = false;
support_map('st_mat_test_prefix') = ['dft_'];

%% No Savings, Vary Borrowing Interest Rate, Single Borrowing Rate (GRID) (Vary fl_r_borr, fl_beta)
%

cl_st_param_keys = {'fl_r_borr', 'fl_beta'};

param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
it_simu_vec_len = 10;
param_tstar_map('fl_r_borr') = linspace(0.01, 0.10, it_simu_vec_len);
param_tstar_map('fl_beta') = linspace(0.90, 0.999, it_simu_vec_len);

% Can not Save
fl_a_max_ori = param_map('fl_a_max');
param_map('fl_a_max') = 0;

% Initialize to be replaced by param_tstar_map inside function
param_map('fl_r_borr') = 0;
param_map('fl_b_bd') = -100;

% Simu Size and Graph Type
it_size_type = 2;
ar_it_plot_sets = [51,52,54, 5,6,153, 61,62,64];
bl_simu_cross = 'g';

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);

close all;

%% Single Saving Interest Rate, Single Borrowing Interest Rate (GRID) (Vary fl_r_borr, fl_beta)


param_map('fl_a_max') = fl_a_max_ori;

cl_st_param_keys = {'fl_r_borr', 'fl_beta'};
param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
it_simu_vec_len = 10;
param_tstar_map('fl_r_borr') = linspace(0.01, 0.10, it_simu_vec_len);
param_tstar_map('fl_beta') = linspace(0.87, 0.97, it_simu_vec_len);

% Simu Size and Graph Type
it_size_type = 2;
ar_it_plot_sets = [51,52,54, 5,6,153, 61,62,64];
bl_simu_cross = 'g';

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);

close all;

%% Single Saving Interest Rate, Single Borrowing Interest Rate (GRID) (Vary fl_r_borr, fl_r_save)

param_map('fl_a_max') = fl_a_max_ori;

cl_st_param_keys = {'fl_r_borr', 'fl_r_save'};
param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
it_simu_vec_len = 10;
param_tstar_map('fl_r_borr') = linspace(0.01, 0.10, it_simu_vec_len);
param_tstar_map('fl_r_save') = linspace(0.00, 0.06, it_simu_vec_len);

% Simu Size and Graph Type
it_size_type = 2;
ar_it_plot_sets = [51,52,54, 5,6,153, 61,62,64];
bl_simu_cross = 'g';

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);

close all;

%% Multiple Borrowing Interest Rates (GRID) (Vary fl_r_borr, fl_r_save)
it_param_set = 9;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% Borrowing Parameters
param_map('bl_default') = bl_default;

% Support Parameters
support_map('bl_replacefile') = false;
support_map('bl_graph_onebyones') = true;
support_map('bl_display_graph_stats') = false;
support_map('st_mat_test_prefix') = ['dft_'];

it_size_type = 10;
ar_it_plot_sets = [51,52,54, 5,6,153, 61,62,64];
bl_simu_cross = 'g';

cl_st_param_keys = {'fl_z_r_borr_poiss_mean', 'fl_r_save'};
it_simu_vec_len = 10;
param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map('fl_z_r_borr_poiss_mean') = linspace(2, 10, it_simu_vec_len);
param_tstar_map('fl_r_save') = linspace(0, 0.06, it_simu_vec_len);

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);

close all;