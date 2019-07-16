%% For + Inf + Borr + Save, Test Shocks *Default* (Save + Borr Distribution)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Given the parameters and bounds here. When default is allowed, at lower
% persistence levels, increasing shock variance does not increase
% consumption variance. But variance for consumption does increase for the
% most persistent shock simulation. Higher persistence also leads to higher
% overall mean levels, more mass on savings. 
%

%% Common Parameters

% Shock Array
ar_fl_z_rho = [0.65, 0.80, 0.95];
ar_fl_z_sig = [0.05, 0.20, 0.35];

% Accuracy
ar_it_a_n = [750];
ar_it_z_n = [15];

% Borrowing/Savings Parameters
bl_default = true;
bl_bridge = true;
bl_rollover = true;
fl_b_bd = -20;
fl_c_min = 0.01; % cmin so low default exists but never chosen
fl_r_fsv = 0.02;
fl_r_inf = 0.065;
fl_r_fbr = 0.045;

% which group to simulate
ar_it_test_grp = [7];
it_simu_vec_len = 4;

%% Set Overriding Parameter Map
param_map = containers.Map('KeyType','char', 'ValueType','any');

% Key Borrowing Controls
param_map('bl_default') = bl_default;
param_map('bl_bridge') = bl_bridge;
param_map('bl_rollover') = bl_rollover;

% Key borrowing float controls
param_map('fl_b_bd') = fl_b_bd;
param_map('fl_c_min') = fl_c_min;

% Interest Rate Simulate Controls
param_map('fl_r_fsv') = fl_r_fsv;
param_map('fl_r_inf') = fl_r_inf;
param_map('fl_r_inf_bridge') = fl_r_inf;
param_map('fl_r_fbr') = fl_r_fbr;

%% Set Support Map

support_map = containers.Map('KeyType','char', 'ValueType','any');
% no overriding changes needed

%% Append to Parameter Arrays
param_test_array_map = containers.Map('KeyType','char', 'ValueType','any');
param_test_array_map('ar_fl_z_sig') = ar_fl_z_sig;

%% Simulate Model with Low Persistence

param_map('fl_beta') = ar_fl_z_rho(1);
fsi_ipwkbzr_fibs_ds_vecsv_testmain(param_map, support_map, ar_it_test_grp)

%% Simulate Model with Medium Persistence

param_map('fl_beta') = ar_fl_z_rho(2);
fsi_ipwkbzr_fibs_ds_vecsv_testmain(param_map, support_map, ar_it_test_grp)

%% Simulate Model with High Persistence

param_map('fl_beta') = ar_fl_z_rho(3);
fsi_ipwkbzr_fibs_ds_vecsv_testmain(param_map, support_map, ar_it_test_grp)
