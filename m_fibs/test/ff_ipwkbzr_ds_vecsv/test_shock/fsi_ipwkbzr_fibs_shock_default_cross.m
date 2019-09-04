%% Test Production Function Parameters (CROSS)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%% Setting of Parameters Locally

close all;
% Set which to graph, simulate over which variables
cl_st_param_keys = {'fl_z_wage_rho', 'fl_z_wage_sig'};

% Generate Arrays of Parameter Values to Loop Over
it_simu_vec_len = 15;
param_tstar_map_local = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map_local('fl_z_wage_rho') = linspace(0, 0.985, it_simu_vec_len);
param_tstar_map_local('fl_z_wage_sig') = linspace(0.01, 0.50, it_simu_vec_len);

%% Set and Load Parameters

bl_default = true;
bl_simu_cross = 'c';

%% Parameter Simulation Arrays

test_map = fsi_ipwkbzr_fibs_ds_support();
ar_it_plot_map = test_map('ar_it_plot_map');
param_tstar_map = test_map('param_tstar_map');
ar_it_size_type = test_map('ar_it_size_type');
ar_it_size_type = [3]
bl_close_all = test_map('bl_close_all');

param_tstar_map = [param_tstar_map; param_tstar_map_local];

%% Parameter Setting Setting

% Generate Benchmark Parameters
it_param_set = 9;
[param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);

% Borrowing Parameters
param_map('bl_default') = bl_default;

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
    ff_az_test_analyze( ...
        ar_it_plot_sets, bl_simu_cross, it_size_type, cl_st_param_keys, ...
        param_map, support_map, param_tstar_map);
    
    if (bl_close_all)
        close all;
    end
    
end