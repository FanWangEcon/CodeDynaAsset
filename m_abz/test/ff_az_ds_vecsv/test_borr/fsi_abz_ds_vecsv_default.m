%% Test Borrowing Default (Save + Borr Distribution)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% @seealso
%
% * test speed: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_speed/html/fsi_az_ds_vecsv_speed.html fsi_az_ds_vecsv_speed>
% * test joint *RANDOM*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_joint/html/fsi_az_ds_vecsv_joint_rand.html fsi_az_ds_vecsv_joint_rand>
% * test interest rate no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc.html fsi_abz_ds_vecsv_nbc>
% * test interest rate no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc_cross.html fsi_abz_ds_vecsv_nbc_cross>
% * test interest rate default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default.html fsi_abz_ds_vecsv_default>
% * test interest rate default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default_cross.html fsi_abz_ds_vecsv_default_cross>
% * test shock default (very low cmin): <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_default_lowcmin.html fsi_abz_ds_vecsv_shk_default_lowcmin>
% * test shock no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_nbc.html fsi_abz_ds_vecsv_shk_nbc>
% * test shock no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_nbc_cross.html fsi_abz_ds_vecsv_shk_nbc_cross>
% * test shock default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_default.html fsi_abz_ds_vecsv_shk_default>
% * test shock default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_default_cross.html fsi_abz_ds_vecsv_shk_default_cross>
% * test preference no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_nbc.html fsi_abz_ds_vecsv_pref_nbc>
% * test preference no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_nbc_cross.html fsi_abz_ds_vecsv_pref_nbc_cross>
% * test preference default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_default.html fsi_abz_ds_vecsv_pref_default>
% * test preference default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_default_cross.html fsi_abz_ds_vecsv_pref_default_cross>
% * test preference default (very low cmin): <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_default_lowcmin.html fsi_abz_ds_vecsv_pref_default_lowcmin>
%

%% Set Shared Parameters

% Accuracy
ar_fl_z_r_borr_n_hg = 3:5:17;
ar_it_a_n_hg = zeros([1,length(ar_fl_z_r_borr_n_hg)]) + 750;

%% Test Minimum Consumption with Defaults, Vary Borrow Bound
% As Borrowing bound relaxes, more borrowing. As borrowing bound is
% relaxed, the fraction of people defaulting is increasing. *pYisMINY*
% captures fraction of households defaulting. Note that this fraction was
% close to almost exactly zero when default was not allowed due to the
% natural borrowing constraint.
%
% Much greater consumption variance, is this calculated correctly?, when
% there is relaxed borrowing constraint. Consumption is not calculated
% correctly for borrowing.
%

close all;

it_size = 4;

% default or not
bl_default = true;
ar_bl_default = zeros([1,it_size]) + bl_default;
% cmin
fl_c_min = 0.025;
ar_fl_c_min = zeros([1,it_size]) + fl_c_min;
% borrow bound
ar_fl_b_bd = linspace(-20, -5, it_size);
% ar_fl_b_bd = zeros([1,it_size]) + fl_b_bd;
% fl_z_r_borr_poiss_mean
fl_z_r_borr_poiss_mean = 0.10;
% ar_fl_z_r_borr_poiss_mean = linspace(0.20, 0.025, it_size);
ar_fl_z_r_borr_poiss_mean = zeros([1,it_size]) + fl_z_r_borr_poiss_mean;

% Simulate Model with Discount = 0.85
for it_cur_param = 1:1:length(ar_bl_default)

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['ar_bl_default = ' num2str(ar_bl_default(it_cur_param))]);
    disp(['ar_fl_c_min = ' num2str(ar_fl_c_min(it_cur_param))]);
    disp(['ar_fl_b_bd = ' num2str(ar_fl_b_bd(it_cur_param))]);
    disp(['ar_fl_z_r_borr_poiss_mean = ' num2str(ar_fl_z_r_borr_poiss_mean(it_cur_param))]);
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
    param_map('fl_z_r_borr_poiss_mean') = ar_fl_z_r_borr_poiss_mean(it_cur_param);

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    support_map('bl_graph_funcgrids') = false;

    for it_accuracy = 1:length(ar_it_a_n_hg)
        % Accuracy Regular
        param_map('it_a_n') = ar_it_a_n_hg(it_accuracy);
        param_map('fl_z_r_borr_n') = ar_fl_z_r_borr_n_hg(it_accuracy);
        it_z_n = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
        param_map('it_z_n') = it_z_n;
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp(['it_a_n = ' num2str(ar_it_a_n_hg(it_accuracy)) ', it_z_n = ' num2str(it_z_n)]);
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
        [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);
        % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
        result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);
        % Call Distribution CProgram
        result_map = ff_az_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
    end


    % Snap
    snapnow;

end

% close all
close all;

%% Test Minimum Consumption with Defaults, Vary Cmin
% A minimum consumption value increases, default increases. As Cmin
% decreases, households still borrow, but there is eventually 0 percent of
% households that are going to the defaulting bound.

close all;

it_size = 4;

% default or not
bl_default = true;
ar_bl_default = zeros([1,it_size]) + bl_default;
% cmin
ar_fl_c_min = linspace(0.1, 0.001, it_size);
% ar_fl_c_min = zeros([1,it_size]) + fl_c_min;
% borrow bound
fl_b_bd = -12.5;
ar_fl_b_bd = zeros([1,it_size]) + fl_b_bd;
% fl_z_r_borr_poiss_mean
fl_z_r_borr_poiss_mean = 0.10;
ar_fl_z_r_borr_poiss_mean = zeros([1,it_size]) + fl_z_r_borr_poiss_mean;

% Simulate Model with Discount = 0.85
for it_cur_param = 1:1:length(ar_bl_default)

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['ar_bl_default = ' num2str(ar_bl_default(it_cur_param))]);
    disp(['ar_fl_c_min = ' num2str(ar_fl_c_min(it_cur_param))]);
    disp(['ar_fl_b_bd = ' num2str(ar_fl_b_bd(it_cur_param))]);
    disp(['ar_fl_z_r_borr_poiss_mean = ' num2str(ar_fl_z_r_borr_poiss_mean(it_cur_param))]);
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
    param_map('fl_z_r_borr_poiss_mean') = ar_fl_z_r_borr_poiss_mean(it_cur_param);

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    support_map('bl_graph_funcgrids') = false;

    for it_accuracy = 1:length(ar_it_a_n_hg)
        % Accuracy Regular
        param_map('it_a_n') = ar_it_a_n_hg(it_accuracy);
        param_map('fl_z_r_borr_n') = ar_fl_z_r_borr_n_hg(it_accuracy);
        it_z_n = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
        param_map('it_z_n') = it_z_n;
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp(['it_a_n = ' num2str(ar_it_a_n_hg(it_accuracy)) ', it_z_n = ' num2str(it_z_n)]);
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
        [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);
        % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
        result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);
        % Call Distribution CProgram
        result_map = ff_az_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
    end

    % Snap
    snapnow;

end

% close all
close all;

%% Test Minimum Consumption with Defaults, Vary Interest rate
% As borrowing interest rate increases, borrowing decreases

close all;

it_size = 10;

% default or not
bl_default = true;
ar_bl_default = zeros([1,it_size]) + bl_default;
% cmin
fl_c_min = 0.035;
ar_fl_c_min = zeros([1,it_size]) + fl_c_min;
% borrow bound
fl_b_bd = -12.5;
ar_fl_b_bd = zeros([1,it_size]) + fl_b_bd;
% fl_z_r_borr_poiss_mean
ar_fl_z_r_borr_poiss_mean = linspace(5, 20, it_size);

% Simulate Model with Discount = 0.85
for it_cur_param = 1:1:length(ar_bl_default)

    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['ar_bl_default = ' num2str(ar_bl_default(it_cur_param))]);
    disp(['ar_fl_c_min = ' num2str(ar_fl_c_min(it_cur_param))]);
    disp(['ar_fl_b_bd = ' num2str(ar_fl_b_bd(it_cur_param))]);
    disp(['ar_fl_z_r_borr_poiss_mean = ' num2str(ar_fl_z_r_borr_poiss_mean(it_cur_param))]);
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
    param_map('fl_z_r_borr_poiss_mean') = ar_fl_z_r_borr_poiss_mean(it_cur_param);

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    support_map('bl_graph_funcgrids') = false;

    for it_accuracy = 1:length(ar_it_a_n_hg)
        % Accuracy Regular
        param_map('it_a_n') = ar_it_a_n_hg(it_accuracy);
        param_map('fl_z_r_borr_n') = ar_fl_z_r_borr_n_hg(it_accuracy);
        it_z_n = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
        param_map('it_z_n') = it_z_n;
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp(['it_a_n = ' num2str(ar_it_a_n_hg(it_accuracy)) ', it_z_n = ' num2str(it_z_n)]);
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        % Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
        [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);
        % Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
        result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);
        % Call Distribution CProgram
        result_map = ff_az_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
    end

    % Snap
    snapnow;

end

% close all
close all;
