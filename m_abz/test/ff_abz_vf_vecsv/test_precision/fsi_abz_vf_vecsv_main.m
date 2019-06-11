%% Tests the ABZ_VF_VECSV Algorithm with several simulations
% 1. benchmark: benchmark not very accurate but useful for comparing against latnerative algorithms
% 2. operational
% 3. most precise

%% Benchmark Simulation, used in Main Testing Files
% benchmark not very accurate but useful for comparing against latnerative algorithms

it_param_set = 4;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% Simulation Accuracy
param_map('it_a_n') = 100;
param_map('it_z_n') = 11;

% Display Parameters
support_map('bl_display') = false;
support_map('bl_display_final') = false;
support_map('bl_time') = true;
support_map('bl_profile') = false;

% Call Program
ff_abz_vf_vecsv(param_map, support_map);


%% Operational Simulation
% fast and accurate enough as main simulation parameters

it_param_set = 4;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% Simulation Accuracy
param_map('it_a_n') = 750;
param_map('it_z_n') = 15;

% Display Parameters
support_map('bl_display') = false;
support_map('bl_display_final') = false;
support_map('bl_time') = true;
support_map('bl_profile') = false;

% Call Program
ff_abz_vf_vecsv(param_map, support_map);

%% High Precision Simulation

it_param_set = 4;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% Simulation Accuracy
param_map('it_a_n') = 2250;
param_map('it_z_n') = 27;

% Display Parameters
support_map('bl_display') = false;
support_map('bl_display_final') = false;
support_map('bl_time') = true;
support_map('bl_profile') = false;

% Call Program
ff_abz_vf_vecsv(param_map, support_map);
