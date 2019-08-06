%% Test Borrowing Default (Save + Borr Distribution)
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

%% Loop over Savings Bounds
clear all;
close all;

ar_fl_a_max = linspace(0,50,2);

for fl_a_max = ar_fl_a_max

    %% Generate Data A: Varying Savings Bounds, Compare Two Borrowing Interest Rates along all AZ points
    % Do not allow for savings, Compare V(a,z;low borrow r) vs V(a,z;high
    % borrow r)


    % Number of Borrowing Interest Rates
    it_size = 2;
    % Do not Allow for Savings
    ar_fl_a_max = zeros([1,it_size]) + fl_a_max;
    % fl_z_r_borr_poiss_mean: see fft_gen_discrete_var
    ar_fl_z_r_borr_max = linspace(0.035, 0.065, it_size);

    % default or not
    bl_default = true;
    ar_bl_default = zeros([1,it_size]) + bl_default;
    % cmin
    fl_c_min = 0.01;
    ar_fl_c_min = zeros([1,it_size]) + fl_c_min;
    % borrow bound
    fl_b_bd = -30;
    ar_fl_b_bd = zeros([1,it_size]) + fl_b_bd;

    % Accuracy
    ar_fl_z_r_borr_n_hg = 1;
    ar_it_a_n_hg = zeros([1,length(ar_fl_z_r_borr_n_hg)]) + 750;
    it_z_wage_n = 15;

    % other parameters
    fl_crra = 1.5;
    fl_beta = 0.96;
    fl_r_save = 0.025;

    % Results Map Collection
    cl_result_map = [];
    it_cl_res_map_ctr = 0;

    % Simulate Model
    for it_cur_param = 1:1:length(ar_bl_default)

        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp(['bl_default = ' num2str(ar_bl_default(it_cur_param))]);
        disp(['fl_c_min = ' num2str(ar_fl_c_min(it_cur_param))]);
        disp(['fl_b_bd = ' num2str(ar_fl_b_bd(it_cur_param))]);
        disp(['fl_a_max = ' num2str(ar_fl_a_max(it_cur_param))]);
        disp(['z_r_borr_max = ' num2str(ar_fl_z_r_borr_max(it_cur_param))]);
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp('');
        disp('');
        disp('');
        disp('');

        % Call Default Parameters <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_set_default_param.html ffs_abz_set_default_param>
        bl_input_override = true;
        it_param_set = 9;
        [param_map, support_map] = ffs_abz_set_default_param(it_param_set);

        % Changing Savings Allowed
        param_map('fl_a_max') = fl_a_max;

        % Borrowing Interest Rates, Same borrow and save rate, no random rate
        param_map('fl_z_r_borr_max') = ar_fl_z_r_borr_max(it_cur_param);
        param_map('fl_z_r_borr_min') = ar_fl_z_r_borr_max(it_cur_param);

        % Borrowing Parameters
        param_map('bl_default') = ar_bl_default(it_cur_param);
        param_map('fl_b_bd') = ar_fl_b_bd(it_cur_param);
        param_map('fl_c_min') = ar_fl_c_min(it_cur_param);

        param_map('fl_beta') = fl_beta;

        % Display Parameters
        support_map('bl_display') = false;
        support_map('bl_display_final') = false;
        support_map('bl_time') = true;
        support_map('bl_profile') = false;

        support_map('bl_graph_funcgrids') = false;

        for it_accuracy = 1:length(ar_it_a_n_hg)
            % Accuracy Regular
            param_map('it_a_n') = ar_it_a_n_hg(it_accuracy);
            param_map('it_z_wage_n') = it_z_wage_n;
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
            % Collect
            it_cl_res_map_ctr = it_cl_res_map_ctr + 1;
            cl_result_map{it_cl_res_map_ctr} = result_map;
        end
        % Snap
        snapnow;

    end

    %% Graph A1: Value Graphs Difference in V(A,Z) for Lower and Higher Borrowing Rate
    % Left to right, increasing z, top to bottom, increasing savings. no
    % borrowing, savings only model, lower borrowing interest rate V(a,z)
    % uniformly higher than higher borrowing interest rate V(a,z).

    close all;

    cl_mt_pol_c = cl_result_map{1}('cl_mt_pol_c');
    cl_mt_pol_c{1};

    mt_val_1 = cl_result_map{1}('mt_val');
    mt_val_2 = cl_result_map{2}('mt_val');

    st_y_title = ['V(a,z;r borrow=' num2str(ar_fl_z_r_borr_max(1)) ') - V(a,z;r borrow=' num2str(ar_fl_z_r_borr_max(2)) ')'];

    figure();
    subplot(1,3,1);
    contour = contourf(mt_val_1, 10);
    clabel(contour);
    xlabel('shock grid points');
    ylabel('asset grid points');
    zlabel('value');
    title(['r borrow=' num2str(ar_fl_z_r_borr_max(1))]);

    subplot(1,3,2);
    contour = contourf(mt_val_2, 10);
    clabel(contour);
    xlabel('shock grid points');
    ylabel('asset grid points');
    zlabel('value');
    title(['r borrow=' num2str(ar_fl_z_r_borr_max(2))]);

    subplot(1,3,3);
    mt_val_diff = mt_val_1 - mt_val_2;

    contour = contourf(mt_val_1 - mt_val_2);
    clabel(contour);
    xlabel('shock grid points');
    ylabel(['asset grid, ' st_y_title]);
    zlabel('value diff');
    title(['Save <= ' num2str(fl_a_max)]);

    %% Graph A2: Distribution Graphs Difference in f(A,Z) for Lower and Higher Borrowing Rate
    % Left to right, increasing z, top to bottom, increasing savings. As
    % interest rate decreases, more mass at lower a points. Households are more
    % likely to borrow. Overall EV is a weighted sum of the previous and this
    % graph.

    mt_dist_1 = cl_result_map{1}('mt_dist');
    mt_dist_2 = cl_result_map{2}('mt_dist');

    st_y_title = ['F(a,z;r borrow=' num2str(ar_fl_z_r_borr_max(1)) ') - F(a,z;r borrow=' num2str(ar_fl_z_r_borr_max(2)) ')'];

    figure();
    subplot(1,3,1);
    surf(mt_dist_1);
    xlabel('shock grid points');
    ylabel('asset grid points');
    zlabel('dist');
    title(['r br=' num2str(ar_fl_z_r_borr_max(1))]);

    subplot(1,3,2);
    surf(mt_dist_2);
    xlabel('shock grid points');
    ylabel('asset grid points');
    zlabel('dist');
    title(['r br=' num2str(ar_fl_z_r_borr_max(2))]);

    subplot(1,3,3);
    contour = contourf(mt_dist_1 - mt_dist_2);
    clabel(contour);
    xlabel('shock grid points');
    ylabel(['asset grid, ' st_y_title]);
    zlabel('dist diff');
    title(['Save <= ' num2str(fl_a_max)]);
end
