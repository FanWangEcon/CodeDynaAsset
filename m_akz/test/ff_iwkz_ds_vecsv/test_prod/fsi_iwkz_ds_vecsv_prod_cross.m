%% Test Production Function (Risky + Safe Asets + Two-Step + Interpolated Distribution) (cross)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.* *cross test*: given (x,y), vary x along X, y along Y,
% one at a time
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
% Note: as for example _fl_beta_ increases, willingness to save increases,
% leading to higher savings, which will exceed the benchmark max grid
% point. So to allow for higher beta, dramatically higher max savings bound
% is needed.

% Set which to graph, simulate over which variables
ar_it_plot_sets = [1,2,101,3,4,102,5,6,103,104,105,106];
bl_simu_cross = 'c';
cl_st_param_keys = {'fl_alpha', 'fl_delta'};

% Generate Benchmark Parameters
it_param_set = 9;
[param_map, support_map] = ffs_akz_set_default_param(it_param_set);
support_map('bl_replacefile') = false;
support_map('bl_graph_onebyones') = true;
support_map('bl_display_graph_stats') = false;

% Generate Arrays of Parameter Values to Loop Over
it_simu_vec_len = 50;
param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map('fl_alpha') = linspace(0.30, 0.50, it_simu_vec_len);
param_tstar_map('fl_delta') = linspace(0.02, 0.14, it_simu_vec_len);

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
