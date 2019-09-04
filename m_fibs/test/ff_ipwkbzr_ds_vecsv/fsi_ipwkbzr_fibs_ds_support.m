%% Shared function for IPWKBZR FISB Testing
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function test_map = fsi_ipwkbzr_fibs_ds_support(varargin)
%% Function with some shared parameters etc. 
% Set commonly shared parameters across testing function here to avoid
% retyping etc. 

%% Parse Parameters

it_simu_vec_len = 15;

default_params = {it_simu_vec_len};
params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};
it_simu_vec_len = default_params{1};

%% What size to Simulate
% 1 is quick, 2 is dense, and 3 is full. s

ar_it_size_type = [1,2,3];

%% Which and What to Plot

ar_it_plot_map = containers.Map('KeyType','char', 'ValueType','any');

ar_it_plot_map('ar_it_plot_sets_cross_limited') = [3,4, 104,106];
ar_it_plot_map('ar_it_plot_sets_cross_full') = [3,4,102,152, 104,106,201,101, 1001,1002,1003,10];

ar_it_plot_map('ar_it_plot_sets_grid_limited') = [51,52,53,54, 5,6,103,153];
ar_it_plot_map('ar_it_plot_sets_grid_full') = [51,52,53,54, 5,6,103,153, 61,62,63,64];

%% Parameter Grids

param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');

param_tstar_map('fl_r_fsv') = linspace(0.00,    0.05, it_simu_vec_len);
param_tstar_map('fl_r_fbr') = linspace(0.03, 0.08, it_simu_vec_len);
param_tstar_map('fl_forbrblk_gap') = linspace(-0.1, -2.5, it_simu_vec_len);

param_tstar_map('fl_z_r_infbr_poiss_mean') = linspace(2, 10, it_simu_vec_len);
param_tstar_map('fl_z_r_infbr_min') = linspace(0.0, 0.025, it_simu_vec_len);
param_tstar_map('fl_z_r_infbr_max') = linspace(0.065, 0.30, it_simu_vec_len);

param_tstar_map('fl_c_min') = linspace(0.001, 0.05, it_simu_vec_len);
param_tstar_map('fl_w') = linspace(0, 0.5, it_simu_vec_len);

param_tstar_map('fl_crra') = linspace(1, 5, it_simu_vec_len);
param_tstar_map('fl_beta') = linspace(0.87, 0.97, it_simu_vec_len);

param_tstar_map('fl_alpha') = linspace(0.30, 0.50, it_simu_vec_len);
param_tstar_map('fl_delta') = linspace(0.02, 0.14, it_simu_vec_len);

%% Parameter Grids requiring Explanations

% When the persistence of shock exceeds 0.9, results investments,
% consumptions all increase dramatically
param_tstar_map('fl_z_wage_rho') = linspace(0, 0.985, it_simu_vec_len);
param_tstar_map('fl_z_wage_sig') = linspace(0.01, 0.50, it_simu_vec_len);


%% Various Controls

bl_close_all = true;

%% All to Test map
test_map = containers.Map('KeyType','char', 'ValueType','any');
test_map('ar_it_plot_map') = ar_it_plot_map;
test_map('param_tstar_map') = param_tstar_map;
test_map('ar_it_size_type') = ar_it_size_type;
test_map('bl_close_all') = bl_close_all;

end
