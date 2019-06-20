%% Test Shock Persistence and Variance (Save + Borr Distribution)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Testing the
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html
% ff_az_ds> program for solving the savings only dynamic
% programming problem.
%
% defaults in ffs_abz_set_default_param.m are:
%
% * param_map('fl_z_rho') = 0.8;
% * param_map('fl_z_sig') = 0.2;
%
% here test three levels of persistence:
%
% * iid shocks
% * 0.5 persistence
% * 0.995 persistence
%
% for each shock, thest at these standard deviations of the log normal
% shock:
%
% * 0.01
% * 0.05
% * 0.30
%
% @seealso
%
% * _PREFERENCE_: savings only quick vs benchmark testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds/test_pref/html/fsi_abz_ds_pref.html fsi_abz_ds_pref>
% * _SHOCK_: savings only asset grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds/test_shock/html/fsi_abz_ds_shock.html fsi_abz_ds_shock>
% * _PREFERENCE_: borr + save quick vs benchmark testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds/test_pref/html/fsi_abz_ds_pref.html fsi_abz_ds_pref>
% * _SHOCK_: borr + save only asset grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds/test_shock/html/fsi_abz_ds_shock.html fsi_abz_ds_shock>
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

    % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_set_default_param.html ffs_abz_set_default_param>
    bl_input_override = true;
    it_param_set = 9;
    [param_map, support_map] = ffs_abz_set_default_param(it_param_set);

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

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
    [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
    result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_az_ds(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

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

    % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_set_default_param.html ffs_abz_set_default_param>
    bl_input_override = true;
    it_param_set = 9;
    [param_map, support_map] = ffs_abz_set_default_param(it_param_set);

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

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
    [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
    result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_az_ds(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;

%% Simulate Model with schok persistence = 0.995 (very persistent)

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

    % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_set_default_param.html ffs_abz_set_default_param>
    bl_input_override = true;
    it_param_set = 9;
    [param_map, support_map] = ffs_abz_set_default_param(it_param_set);

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

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
    [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
    result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_az_ds(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;
clear all;
