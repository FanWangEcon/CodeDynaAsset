%% Test Borrowing (Save + Borr Distribution)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Testing the
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vec.html
% ff_az_ds_vec> program for solving the savings only dynamic
% programming problem.
%
% No default, change borrowing interest rates. Default, changin minimum
% consumption etc. 
%


%% Test Borrowing Interest Rates, No Default, NBC
% The natural borrowing constraint restricts borrowing. Note that
% *pYisMINY* in tables below is generally very close to 0. That is the
% correct mass. As we shift the formal borrowing interest rate, the natural
% borrowing constraints changes (lower borrowing interest rate, more
% borrowing allowed), but the mass at the borrowing bound is 0 because at
% the bound, households would go to negative infinite utility, approximated
% by -9999 in default parameters pYisMINY. Note the *min* value in table
% shows how the natural borrowing constraint is adjusting. that is the NBC
% bound.
%

close all;
clear all;

it_size = 5;

% default or not
bl_default = false;
ar_bl_default = [false, false, true];
ar_bl_default = zeros([1,it_size]) + bl_default;
% cmin
fl_c_min = 0.001;
ar_fl_c_min = [0.001, 0.005, 0.010];
ar_fl_c_min = linspace(0.001, 0.005, it_size);
ar_fl_c_min = zeros([1,it_size]) + fl_c_min;
% borrow bound
fl_b_bd = -20;
ar_fl_b_bd = [-20, -15, -10];
ar_fl_b_bd = linspace(-20, 0, it_size);
ar_fl_b_bd = zeros([1,it_size]) + fl_b_bd;
% fl_r_borr
fl_r_borr = 0.035;
ar_fl_r_borr = [-999, -999, -999];
ar_fl_r_borr = linspace(0.20, 0.025, it_size);
% ar_fl_r_borr = zeros([1,it_size]) + fl_r_borr;

% other parameters
it_a_n = 750;
it_z_n = 15;
fl_crra = 1.5;
fl_beta = 0.94;
fl_r_save = 0.025;

% Simulate Model

for it_cur_param = 1:1:length(ar_bl_default)

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['ar_bl_default = ' num2str(ar_bl_default(it_cur_param))]);
    disp(['ar_fl_c_min = ' num2str(ar_fl_c_min(it_cur_param))]);
    disp(['ar_fl_b_bd = ' num2str(ar_fl_b_bd(it_cur_param))]);
    disp(['ar_fl_r_borr = ' num2str(ar_fl_r_borr(it_cur_param))]);
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

    % Borrowing Parameters
    param_map('bl_default') = ar_bl_default(it_cur_param);
    param_map('fl_c_min') = ar_fl_c_min(it_cur_param);
    param_map('fl_b_bd') = ar_fl_b_bd(it_cur_param);
    
    % Interest Rates
    param_map('fl_r_save') = fl_r_save;
    param_map('fl_r_borr') = ar_fl_r_borr(it_cur_param);
    
    % Simulation Accuracy
    param_map('it_a_n') = it_a_n;
    param_map('it_z_n') = it_z_n;
    param_map('fl_crra') = fl_crra;
    param_map('fl_beta') = fl_beta;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;
    
    support_map('bl_graph_funcgrids') = false;

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
    [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
    result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_az_ds_vec(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;

%% Test Borrowing Bounds with Defaults
% Borrowing could now go to the borrowing bound, which is the point of
% default. Note the fraction of households at *pYisMINY*, unlike in the
% exercise above is no longer 0, but potentially a positive mass. The
% households at these points are defaulting. A crucial parameter that
% controls default is _fl_c_min_, what is the u(c) utility today when
% households defualt? If this is very low, then even if default is allowed,
% it will not be an option taken by anyone, and *pYisMINY* will be zero.
% Note here that changing the borrowing bound is also impactful. When the
% borrowing bound is tight, with default allowed, households are unwilling
% to default, this is because if you are forced to get u(c_min) having not
% borrowed much, so borrowing is not worth-while.

close all;
clear all;

it_size = 3;

% default or not
bl_default = true;
ar_bl_default = [false, false, true];
ar_bl_default = zeros([1,it_size]) + bl_default;
% cmin
fl_c_min = 0.025;
ar_fl_c_min = [0.001, 0.005, 0.010];
ar_fl_c_min = linspace(0.001, 0.005, it_size);
ar_fl_c_min = zeros([1,it_size]) + fl_c_min;
% borrow bound
fl_b_bd = -20;
ar_fl_b_bd = [-20, -15, -10];
ar_fl_b_bd = linspace(-20, -5, it_size);
% ar_fl_b_bd = zeros([1,it_size]) + fl_b_bd;
% fl_r_borr
fl_r_borr = 0.10;
ar_fl_r_borr = [-999, -999, -999];
% ar_cur_param_fur = linspace(0.20, 0.025, it_size);
ar_fl_r_borr = zeros([1,it_size]) + fl_r_borr;

% other parameters
it_a_n = 750;
it_z_n = 15;
fl_crra = 1.5;
fl_beta = 0.94;
fl_r_save = 0.025;

% Simulate Model with Discount = 0.85

for it_cur_param = 1:1:length(ar_bl_default)

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['ar_bl_default = ' num2str(ar_bl_default(it_cur_param))]);
    disp(['ar_fl_c_min = ' num2str(ar_fl_c_min(it_cur_param))]);
    disp(['ar_fl_b_bd = ' num2str(ar_fl_b_bd(it_cur_param))]);
    disp(['ar_fl_r_borr = ' num2str(ar_fl_r_borr(it_cur_param))]);
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

    % Borrowing Parameters
    param_map('bl_default') = ar_bl_default(it_cur_param);
    param_map('fl_c_min') = ar_fl_c_min(it_cur_param);
    param_map('fl_b_bd') = ar_fl_b_bd(it_cur_param);
    
    % Interest Rates
    param_map('fl_r_save') = fl_r_save;
    param_map('fl_r_borr') = ar_fl_r_borr(it_cur_param);
    
    % Simulation Accuracy
    param_map('it_a_n') = it_a_n;
    param_map('it_z_n') = it_z_n;
    param_map('fl_crra') = fl_crra;
    param_map('fl_beta') = fl_beta;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;
    
    support_map('bl_graph_funcgrids') = false;

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
    [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
    result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_az_ds_vec(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;

%% Test Minimum Consumption with Defaults
% A minimum consumption value increases, default increases. 

close all;
clear all;

it_size = 4;

% default or not
bl_default = true;
ar_bl_default = [false, false, true];
ar_bl_default = zeros([1,it_size]) + bl_default;
% cmin
fl_c_min = 0.025;
ar_fl_c_min = [0.001, 0.005, 0.010];
ar_fl_c_min = linspace(0.1, 0.001, it_size);
% ar_fl_c_min = zeros([1,it_size]) + fl_c_min;
% borrow bound
fl_b_bd = -12.5;
ar_fl_b_bd = [-20, -15, -10];
% ar_fl_b_bd = linspace(-20, -5, it_size);
ar_fl_b_bd = zeros([1,it_size]) + fl_b_bd;
% fl_r_borr
fl_r_borr = 0.10;
ar_fl_r_borr = [-999, -999, -999];
% ar_cur_param_fur = linspace(0.20, 0.025, it_size);
ar_fl_r_borr = zeros([1,it_size]) + fl_r_borr;

% other parameters
it_a_n = 750;
it_z_n = 15;
fl_crra = 1.5;
fl_beta = 0.94;
fl_r_save = 0.025;

% Simulate Model with Discount = 0.85

for it_cur_param = 1:1:length(ar_bl_default)

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['ar_bl_default = ' num2str(ar_bl_default(it_cur_param))]);
    disp(['ar_fl_c_min = ' num2str(ar_fl_c_min(it_cur_param))]);
    disp(['ar_fl_b_bd = ' num2str(ar_fl_b_bd(it_cur_param))]);
    disp(['ar_fl_r_borr = ' num2str(ar_fl_r_borr(it_cur_param))]);
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

    % Borrowing Parameters
    param_map('bl_default') = ar_bl_default(it_cur_param);
    param_map('fl_c_min') = ar_fl_c_min(it_cur_param);
    param_map('fl_b_bd') = ar_fl_b_bd(it_cur_param);
    
    % Interest Rates
    param_map('fl_r_save') = fl_r_save;
    param_map('fl_r_borr') = ar_fl_r_borr(it_cur_param);
    
    % Simulation Accuracy
    param_map('it_a_n') = it_a_n;
    param_map('it_z_n') = it_z_n;
    param_map('fl_crra') = fl_crra;
    param_map('fl_beta') = fl_beta;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;
    
    support_map('bl_graph_funcgrids') = false;

    % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
    [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);

    % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
    result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);

    % Call Distribution CProgram
    result_map = ff_az_ds_vec(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Snap
    snapnow;

end

% close all
close all;