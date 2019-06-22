%% Test Preference No Default (Save + Borr Distribution)
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
% * param_map('fl_beta') = 0.94;
% * param_map('fl_crra') = 1.5;
%
% here test three levels of discount:
%
% * 0.94
% * 0.96
% * 0.98
%
% for each shock, thest at these crra levels
%
% * log (1)
% * 1.5
% * 2.0
%
% As crra level increases, general consumption variance goes down, this is
% achieved by increasing asset variance. Mean savings and mean consumption
% both go up. The percentiles of consumption is tighter.
%
% Borrowing is bound by the natural borrowing constraint, which is never
% binding. Zero mass at asset min which is NBC bound. 
%
% As discount factor increases, savings increas and borrowing decrease. At
% 0.98, only the 0.1 percentile of households are participating in
% borrowing. 
%
% With higher discount factor, model solution takes longer to converge.
% Also the case when crra parameter increases. Max speed 7.9 seconds, crra
% = 2 and beta = .98; min speed 1.9 seconds, crra = 1, beta = 0.94
%

%% Set Shared Parameters

close all;
clear all;

ar_fl_beta = [0.94, 0.96, 0.98];
ar_fl_crra = [1, 1.5, 2.0];
it_a_n = 750;
it_z_n = 15;

bl_default = false;
fl_c_min = 0.035; % irrelevant when bl_default = false
fl_b_bd = -20;   
fl_r_save = 0.02;
fl_r_borr = 0.065;

%% Simulate Model with Discount = 0.94

for fl_crra = ar_fl_crra

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_beta = ' num2str(ar_fl_beta(1))]);
    disp(['fl_crra = ' num2str(fl_crra)]);
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
    param_map('fl_beta') = ar_fl_beta(1);
    param_map('fl_crra') = fl_crra;

    % Borrowing Parameters
    param_map('bl_default') = bl_default;
    param_map('fl_c_min') = fl_c_min;
    param_map('fl_b_bd') = fl_b_bd;
    
    % Interest Rates
    param_map('fl_r_save') = fl_r_save;
    param_map('fl_r_borr') = fl_r_borr;
        
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
    result_map = ff_az_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;

%% Simulate Model with Discount = 0.96

close all

for fl_crra = ar_fl_crra

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_beta = ' num2str(ar_fl_beta(2))]);
    disp(['fl_crra = ' num2str(fl_crra)]);
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
    param_map('fl_beta') = ar_fl_beta(2);
    param_map('fl_crra') = fl_crra;

    % Borrowing Parameters
    param_map('bl_default') = bl_default;
    param_map('fl_c_min') = fl_c_min;
    param_map('fl_b_bd') = fl_b_bd;
    
    % Interest Rates
    param_map('fl_r_save') = fl_r_save;
    param_map('fl_r_borr') = fl_r_borr;
    
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

%% Simulate Model with Discount = 0.98

close all

for fl_crra = ar_fl_crra

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_beta = ' num2str(ar_fl_beta(3))]);
    disp(['fl_crra = ' num2str(fl_crra)]);
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
    param_map('fl_beta') = ar_fl_beta(3);
    param_map('fl_crra') = fl_crra;

    % Borrowing Parameters
    param_map('bl_default') = bl_default;
    param_map('fl_c_min') = fl_c_min;
    param_map('fl_b_bd') = fl_b_bd;
    
    % Interest Rates
    param_map('fl_r_save') = fl_r_save;
    param_map('fl_r_borr') = fl_r_borr;
    
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
