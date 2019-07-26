%% Test Shock Persistence and Variance (Risky + Safe Asets + Two-Step + Interpolated Distribution)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_test_analyze.html ff_az_test_analyze>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_akz_set_default_param.html ffs_akz_set_default_param>
%
% @seealso
%
% * _SPEED_ risky + safe (two-step interp) overall benchmark speed testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_speed/html/fsi_ikwz_ds_vecsv_speed.html fsi_ikwz_ds_vecsv_speed>
% * _PREFERENCE_ risky + safe (two-step interp) preference testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_pref/html/fsi_ikwz_ds_vecsv_pref.html fsi_ikwz_ds_vecsv_pref>
% * _PREFERENCE_ risky + safe (two-step interp) preference testing *cross*:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_pref/html/fsi_ikwz_ds_vecsv_pref_cross.html fsi_ikwz_ds_vecsv_pref_cross>
% * _SHOCK_ risky + safe (two-step interp) shock testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_shock/html/fsi_ikwz_ds_vecsv_shock.html fsi_ikwz_ds_vecsv_shock>
% * _SHOCK_ risky + safe (two-step interp) shock testing *cross*:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_shock/html/fsi_ikwz_ds_vecsv_shock_cross.html fsi_ikwz_ds_vecsv_shock_cross>
% * _CAPITAL_ risky + safe (two-step interp) capital return testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_prod/html/fsi_ikwz_ds_vecsv_prod.html fsi_ikwz_ds_vecsv_prod>
% * _CAPITAL_ risky + safe (two-step interp) capital return testing *cross*:  <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_prod/html/fsi_ikwz_ds_vecsv_prod_cross.html fsi_ikwz_ds_vecsv_prod_cross>
% * _PRICE_ risky + safe (two-step interp) wage and interest rate testing *cross*: adjust wage and savings rate
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_price/html/fsi_ikwz_ds_vecsv_price_cross.html fsi_ikwz_ds_vecsv_price_cross>
% * _JOINT_ all parameters random draws joint test
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_joint/html/fsi_ikwz_ds_vecsv_joint_rand.html fsi_ikwz_ds_vecsv_joint_rand>
%

%% Set Shared Parameters

close all;
clear all;

ar_fl_z_rho = [0.0, 0.50, 0.99];
ar_fl_z_sig = [0.05, 0.10, 0.3];
it_a_n = 750;
it_z_n = 15;

%% Simulate Model with schok persistence = 0.0, IID

for fl_z_sig = ar_fl_z_sig

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_z_rho = ' num2str(ar_fl_z_rho(1))]);
    disp(['fl_z_sig = ' num2str(fl_z_sig)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_akz_set_default_param.html ffs_akz_set_default_param>
    bl_input_override = true;
    it_param_set = 9;
    [param_map, support_map] = ffs_akz_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_a_n') = it_a_n;
    param_map('it_z_n') = it_z_n;
    param_map('fl_z_rho') = ar_fl_z_rho(1);
    param_map('fl_z_sig') = fl_z_sig;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;
    support_map('bl_graph_coh_t_coh') = true;

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_akz_get_funcgrid.html ffs_akz_get_funcgrid>
    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_iwkz_vf_vecsv.html ff_iwkz_vf_vecsv>
    result_map = ff_iwkz_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_iwkz_ds_vec(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;

%% Simulate Model with schok persistence = 0.5

close all

for fl_z_sig = ar_fl_z_sig

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_z_rho = ' num2str(ar_fl_z_rho(2))]);
    disp(['fl_z_sig = ' num2str(fl_z_sig)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_akz_set_default_param.html ffs_akz_set_default_param>
    bl_input_override = true;
    it_param_set = 9;
    [param_map, support_map] = ffs_akz_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_a_n') = it_a_n;
    param_map('it_z_n') = it_z_n;
    param_map('fl_z_rho') = ar_fl_z_rho(2);
    param_map('fl_z_sig') = fl_z_sig;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;
    support_map('bl_graph_coh_t_coh') = true;

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_akz_get_funcgrid.html ffs_akz_get_funcgrid>
    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_iwkz_vf_vecsv.html ff_iwkz_vf_vecsv>
    result_map = ff_iwkz_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_iwkz_ds_vec(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;

%% Simulate Model with schok persistence = 0.99 (very persistent)

close all

for fl_z_sig = ar_fl_z_sig

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_z_rho = ' num2str(ar_fl_z_rho(3))]);
    disp(['fl_z_sig = ' num2str(fl_z_sig)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_akz_set_default_param.html ffs_akz_set_default_param>
    bl_input_override = true;
    it_param_set = 9;
    [param_map, support_map] = ffs_akz_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_a_n') = it_a_n;
    param_map('it_z_n') = it_z_n;
    param_map('fl_z_rho') = ar_fl_z_rho(3);
    param_map('fl_z_sig') = fl_z_sig;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;
    support_map('bl_graph_coh_t_coh') = true;

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_akz_get_funcgrid.html ffs_akz_get_funcgrid>
    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_iwkz_vf_vecsv.html ff_iwkz_vf_vecsv>
    result_map = ff_iwkz_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_iwkz_ds_vec(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;
clear all;
