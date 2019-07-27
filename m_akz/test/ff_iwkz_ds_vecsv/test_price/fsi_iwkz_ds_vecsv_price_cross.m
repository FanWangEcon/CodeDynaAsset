%% Test Wage and Savings Rate (Risky + Safe Asets + Two-Step + Interpolated Distribution)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% wage refers to a multiple of minimum income. exogenous wage, unrelated to
% production function.
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_test_analyze.html ff_az_test_analyze>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_default_param.html ffs_akz_set_default_param>
%
% @seealso
%
% * _SPEED_ risky + safe (two-step interp) overall benchmark speed testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_speed/html/fsi_iwkz_ds_vecsv_speed.html fsi_iwkz_ds_vecsv_speed>
% * _PREFERENCE_ risky + safe (two-step interp) preference testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_pref/html/fsi_iwkz_ds_vecsv_pref.html fsi_iwkz_ds_vecsv_pref>
% * _PREFERENCE_ risky + safe (two-step interp) preference testing *cross*:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_pref/html/fsi_iwkz_ds_vecsv_pref_cross.html fsi_iwkz_ds_vecsv_pref_cross>
% * _SHOCK_ risky + safe (two-step interp) shock testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_shock/html/fsi_iwkz_ds_vecsv_shock.html fsi_iwkz_ds_vecsv_shock>
% * _SHOCK_ risky + safe (two-step interp) shock testing *cross*:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_shock/html/fsi_iwkz_ds_vecsv_shock_cross.html fsi_iwkz_ds_vecsv_shock_cross>
% * _CAPITAL_ risky + safe (two-step interp) capital return testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_prod/html/fsi_iwkz_ds_vecsv_prod.html fsi_iwkz_ds_vecsv_prod>
% * _CAPITAL_ risky + safe (two-step interp) capital return testing *cross*:  <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_prod/html/fsi_iwkz_ds_vecsv_prod_cross.html fsi_iwkz_ds_vecsv_prod_cross>
% * _PRICE_ risky + safe (two-step interp) wage and interest rate testing *cross*: adjust wage and savings rate
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_price/html/fsi_iwkz_ds_vecsv_price_cross.html fsi_iwkz_ds_vecsv_price_cross>
% * _JOINT_ all parameters random draws joint test
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_joint/html/fsi_iwkz_ds_vecsv_joint_rand.html fsi_iwkz_ds_vecsv_joint_rand>
%

%% Simulate and Graph
% Note that a shift in wage just rescales the model

% Set which to graph, simulate over which variables
ar_it_plot_sets = [1,2,101,3,4,102,5,6,103,104,105,106];
bl_simu_cross = 'c';
cl_st_param_keys = {'fl_w', 'fl_r_save'};

% Generate Benchmark Parameters
it_param_set = 9;
[param_map, support_map] = ffs_akz_set_default_param(it_param_set);
support_map('bl_replacefile') = false;
support_map('bl_graph_onebyones') = true;
support_map('bl_display_graph_stats') = false;

% Generate Arrays of Parameter Values to Loop Over
it_simu_vec_len = 50;
param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map('fl_w') = linspace(0.25, 5, it_simu_vec_len);
param_tstar_map('fl_r_save') = linspace(0, 0.06, it_simu_vec_len);

%% Medium Grid Simulation
it_size_type = 2;

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);

%% Larger Grid Simulation
it_size_type = 3;

% Simulate along parameters
ff_az_test_analyze( ...
    ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
    param_map, support_map, param_tstar_map);
