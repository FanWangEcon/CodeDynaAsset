%% Test Joint Randomly (Risky + Safe Asset + Save + Borr + R Shock + Interpolated-Percentage), Test Preference, Price and Shocks
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.* *cross test*: given (x,y), vary x along X, y along Y,
% one at a time
%

%% Setting of Parameters Locally

close all;
% Set which to graph, simulate over which variables
cl_st_param_keys = {'fl_r_fsv','fl_r_fbr','fl_forbrblk_gap', ...
                    'fl_z_r_infbr_poiss_mean','fl_z_r_infbr_max', 'fl_z_r_infbr_min',...
                    'fl_c_min','fl_w',...
                    'fl_crra','fl_beta',...
                    'fl_alpha','fl_delta',...
                    'fl_z_wage_rho','fl_z_wage_sig'};

% Generate Arrays of Parameter Values to Loop Over
it_simu_vec_len = 2;
param_tstar_map_local = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map_local('fl_r_fsv') = linspace(0.00, 0.02, it_simu_vec_len);
param_tstar_map_local('fl_r_fbr') = linspace(0.05, 0.09, it_simu_vec_len);
param_tstar_map_local('fl_forbrblk_gap') = linspace(-0.1, -2.5, it_simu_vec_len);

param_tstar_map_local('fl_z_r_infbr_poiss_mean') = linspace(1, 20, it_simu_vec_len);
param_tstar_map_local('fl_z_r_infbr_max') = linspace(0.15, 0.30, it_simu_vec_len);
param_tstar_map_local('fl_z_r_infbr_min') = linspace(0.0, 0.025, it_simu_vec_len);

param_tstar_map_local('fl_c_min') = linspace(0.001, 0.05, it_simu_vec_len);
param_tstar_map_local('fl_w') = linspace(0, 0.5, it_simu_vec_len);

param_tstar_map_local('fl_crra') = linspace(1, 2.5, it_simu_vec_len);
param_tstar_map_local('fl_beta') = linspace(0.87, 0.97, it_simu_vec_len);

param_tstar_map_local('fl_alpha') = linspace(0.10, 0.50, it_simu_vec_len);
param_tstar_map_local('fl_delta') = linspace(0.02, 0.14, it_simu_vec_len);

param_tstar_map_local('fl_z_wage_rho') = linspace(0, 0.985, it_simu_vec_len);
param_tstar_map_local('fl_z_wage_sig') = linspace(0.01, 0.50, it_simu_vec_len);

%% Set and Load Parameters

bl_default = true;
bl_simu_cross = 'r';

%% Parameter Simulation Arrays

test_map = fsi_ipwkbzr_fibs_ds_support();
ar_it_plot_map = test_map('ar_it_plot_map');
param_tstar_map = test_map('param_tstar_map');
ar_it_size_type = test_map('ar_it_size_type');
bl_close_all = test_map('bl_close_all');

param_tstar_map = [param_tstar_map; param_tstar_map_local];

%% Parameter Setting Setting

% Generate Benchmark Parameters
it_param_set = 9;
[param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);

% Borrowing Parameters
param_map('bl_default') = bl_default;
param_map('it_st_simu_type_g_seed') = 123;

% Support Parameters
support_map('bl_replacefile') = false;
support_map('bl_graph_onebyones') = true;
support_map('bl_display_graph_stats') = false;
support_map('st_mat_test_prefix') = ['dft_'];

%% Quick CROSS Simulation (Limited Graphs)

if (ismember(1, ar_it_size_type))

    it_size_type = 1;
    ar_it_plot_sets = ar_it_plot_map('ar_it_plot_sets_cross_limited');

    % Simulate along parameters
    param_map('it_st_simu_type_g_simun') = 4000;
    ff_az_test_analyze( ...
        ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
        param_map, support_map, param_tstar_map);

    if (bl_close_all)
        close all;
    end

end

%% Medium CROSS Simulation (Limited Graphs)

if (ismember(2, ar_it_size_type))

    it_size_type = 2;
    ar_it_plot_sets = ar_it_plot_map('ar_it_plot_sets_cross_limited');

    % Simulate along parameters
    param_map('it_st_simu_type_g_simun') = 1000;
    ff_az_test_analyze( ...
        ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
        param_map, support_map, param_tstar_map);

    if (bl_close_all)
        close all;
    end

end

%% Large CROSS Simulation

if (ismember(3, ar_it_size_type))

    it_size_type = 3;
    ar_it_plot_sets = ar_it_plot_map('ar_it_plot_sets_cross_full');

    % Simulate along parameters
    param_map('it_st_simu_type_g_simun') = 250;    
    ff_az_test_analyze( ...
        ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
        param_map, support_map, param_tstar_map);

    if (bl_close_all)
        close all;
    end

end
