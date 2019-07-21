%% X and Y array simulations fixing x and y
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%% Simulate Along Preference and Price
% specify which parameters to be cross-simulated over, and the parameter
% arrays. 

[param_map, support_map] = ffs_az_set_default_param(it_param_set);

bl_simu_cross = true;
cl_st_param_keys = {'fl_crra', 'fl_beta', 'fl_w', 'fl_r_save'};
it_simu_vec_len = 5;
it_size_type = 1;

param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map('fl_crra') = {'crra', linspace(1, 2, it_simu_vec_len)};
param_tstar_map('fl_beta') = {'discount', linspace(0.94, 0.98, it_simu_vec_len)};
param_tstar_map('fl_w') = {'wage' , linspace(1.1, 1.4, it_simu_vec_len)};
param_tstar_map('fl_r_save') = {'save rate', linspace(0.01, 0.04, it_simu_vec_len)};

ff_az_test_gen(bl_simu_cross, it_size_type, cl_st_param_keys, param_map, support_map, param_tstar_map);
