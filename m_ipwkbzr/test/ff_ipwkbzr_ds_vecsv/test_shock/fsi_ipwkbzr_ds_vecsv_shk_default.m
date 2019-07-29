%% Test Shock Persistence and Variance *Default* (Risky + Safe Asset + Save + Borr + R Shock + Interpolated-Percentage)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Compare against
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_iwkz_ds_vecsv/test_shock/html/fsi_iwkz_ds_vecsv_shock.html
% fsi_iwkz_ds_vecsv_shock> from *iwkz*. Here, we use *ipwkz*.
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_test_analyze.html ff_az_test_analyze>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/paramfunc/html/ffs_ipwkbzr_set_default_param.html ffs_ipwkbzr_set_default_param>
%
% @seealso
%
% * test speed: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_speed/html/fsi_ipwkbzr_ds_vecsv_speed_default.html fsi_ipwkbzr_ds_vecsv_speed_default>
% * test joint *RANDOM*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_joint/html/fsi_ipwkbzr_ds_vecsv_joint_default_rand.html fsi_ipwkbzr_ds_vecsv_joint_default_rand>
%
% * test interest rate default: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_borr/html/fsi_ipwkbzr_ds_vecsv_default.html fsi_ipwkbzr_ds_vecsv_default>
% * test interest rate default *RANDOM*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_borr/html/fsi_ipwkbzr_ds_vecsv_default_rand.html fsi_ipwkbzr_ds_vecsv_default_rand>
% * test interest rate default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_borr/html/fsi_ipwkbzr_ds_vecsv_default_cross.html fsi_ipwkbzr_ds_vecsv_default_cross>
% * test interest rate no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_borr/html/fsi_ipwkbzr_ds_vecsv_nbc_cross.html fsi_ipwkbzr_ds_vecsv_nbc_cross>
%
% * test preference default: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_pref/html/fsi_ipwkbzr_ds_vecsv_pref_default.html fsi_ipwkbzr_ds_vecsv_pref_default>
% * test preference default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_pref/html/fsi_ipwkbzr_ds_vecsv_pref_default_cross.html fsi_ipwkbzr_ds_vecsv_pref_default_cross>
% * test preference no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_pref/html/fsi_ipwkbzr_ds_vecsv_pref_nbc_cross.html fsi_ipwkbzr_ds_vecsv_pref_nbc_cross>
%
% * test ymin and save r default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_price/html/fsi_ipwkbzr_ds_vecsv_price_default_cross.html fsi_ipwkbzr_ds_vecsv_price_default_cross>
% * test ymin and save r no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_price/html/fsi_ipwkbzr_ds_vecsv_price_nbc_cross.html fsi_ipwkbzr_ds_vecsv_price_nbc_cross>
%
% * test preference default: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_prod/html/fsi_ipwkbzr_ds_vecsv_prod_default.html fsi_ipwkbzr_ds_vecsv_prod_default>
% * test production default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_prod/html/fsi_ipwkbzr_ds_vecsv_prod_default_cross.html fsi_ipwkbzr_ds_vecsv_prod_default_cross>
% * test production no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_prod/html/fsi_ipwkbzr_ds_vecsv_prod_nbc_cross.html fsi_ipwkbzr_ds_vecsv_prod_nbc_cross>
%
% * test shock default: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_shock/html/fsi_ipwkbzr_ds_vecsv_shk_default.html fsi_ipwkbzr_ds_vecsv_shk_default>
% * test shock default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_shock/html/fsi_ipwkbzr_ds_vecsv_shk_default_cross.html fsi_ipwkbzr_ds_vecsv_shk_default_cross>
% * test shock no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/test/ff_ipwkbzr_ds_vecsv/test_shock/html/fsi_ipwkbzr_ds_vecsv_shk_nbc_cross.html fsi_ipwkbzr_ds_vecsv_shk_nbc_cross>
%

%% Set Shared Parameters

close all;
clear all;

ar_fl_z_wage_rho = [0.0, 0.50, 0.99];
ar_fl_z_wage_sig = [0.05, 0.10, 0.3];

%% Simulate Model with schok persistence = 0.0, IID

for fl_z_wage_sig = ar_fl_z_wage_sig

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_z_wage_rho = ' num2str(ar_fl_z_wage_rho(1))]);
    disp(['fl_z_wage_sig = ' num2str(fl_z_wage_sig)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/paramfunc/html/ffs_ipwkbzr_set_default_param.html ffs_ipwkbzr_set_default_param>
    bl_input_override = true;
    it_param_set = 9;
    [param_map, support_map] = ffs_ipwkbzr_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('fl_z_wage_rho') = ar_fl_z_wage_rho(1);
    param_map('fl_z_wage_sig') = fl_z_wage_sig;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;
    support_map('bl_graph_coh_t_coh') = true;

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_ipwkbzr_get_funcgrid.html ffs_ipwkbzr_get_funcgrid>
    [armt_map, func_map] = ffs_ipwkbzr_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_ipwkbzr_vf_vecsv.html ff_ipwkbzr_vf_vecsv>
    result_map = ff_ipwkbzr_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;

%% Simulate Model with schok persistence = 0.5

close all

for fl_z_wage_sig = ar_fl_z_wage_sig

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_z_wage_rho = ' num2str(ar_fl_z_wage_rho(2))]);
    disp(['fl_z_wage_sig = ' num2str(fl_z_wage_sig)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/paramfunc/html/ffs_ipwkbzr_set_default_param.html ffs_ipwkbzr_set_default_param>
    bl_input_override = true;
    it_param_set = 9;
    [param_map, support_map] = ffs_ipwkbzr_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('fl_z_wage_rho') = ar_fl_z_wage_rho(2);
    param_map('fl_z_wage_sig') = fl_z_wage_sig;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;
    support_map('bl_graph_coh_t_coh') = true;

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_ipwkbzr_get_funcgrid.html ffs_ipwkbzr_get_funcgrid>
    [armt_map, func_map] = ffs_ipwkbzr_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_ipwkbzr_vf_vecsv.html ff_ipwkbzr_vf_vecsv>
    result_map = ff_ipwkbzr_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;

%% Simulate Model with schok persistence = 0.99 (very persistent)

close all

for fl_z_wage_sig = ar_fl_z_wage_sig

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_z_wage_rho = ' num2str(ar_fl_z_wage_rho(3))]);
    disp(['fl_z_wage_sig = ' num2str(fl_z_wage_sig)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbzr/paramfunc/html/ffs_ipwkbzr_set_default_param.html ffs_ipwkbzr_set_default_param>
    bl_input_override = true;
    it_param_set = 9;
    [param_map, support_map] = ffs_ipwkbzr_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('fl_z_wage_rho') = ar_fl_z_wage_rho(3);
    param_map('fl_z_wage_sig') = fl_z_wage_sig;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;
    support_map('bl_graph_coh_t_coh') = true;

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_ipwkbzr_get_funcgrid.html ffs_ipwkbzr_get_funcgrid>
    [armt_map, func_map] = ffs_ipwkbzr_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_ipwkbzr_vf_vecsv.html ff_ipwkbzr_vf_vecsv>
    result_map = ff_ipwkbzr_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;
clear all;
