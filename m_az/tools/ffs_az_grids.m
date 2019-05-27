function [armt_map] = ffs_az_grids(varargin)
% Obtain Matrixes/Vectors of Choices

%% Catch Error
params_len = length(varargin);
if params_len > 2
    error('grid_lambda:TooManyOptionalParameters', ...
        'allows at most 2 optional parameters');
end

%% Defaults
[param_map, support_map] = ffs_az_default_param();
support_map('bl_graph_grids') = true;
support_map('bl_display_grids') = true;

default_maps = {param_map, support_map};

%% Parse Parameters
% numvarargs is the number of varagin inputted
[default_maps{1:params_len}] = varargin{:};
param_map = [param_map; default_maps{1}];
support_map = [support_map; default_maps{2}];

params_group = values(param_map, {'it_z_n', 'fl_z_mu', 'fl_z_rho', 'fl_z_sig'});
[it_z_n, fl_z_mu, fl_z_rho, fl_z_sig] = params_group{:};

params_group = values(param_map, {'fl_b_bd', 'fl_a_min', 'fl_a_max', 'bl_loglin', 'fl_loglin_threshold', 'it_a_n'});
[fl_b_bd, fl_a_min, fl_a_max, bl_loglin, fl_loglin_threshold, it_a_n] = params_group{:};

params_group = values(support_map, {'bl_graph_grids', 'bl_display_grids'});
[bl_graph_grids, bl_display_grids] = params_group{:};

%% Asset and Choice Grid
if (bl_loglin)
    % C:\Users\fan\M4Econ\asset\grid\ff_grid_loglin.m
    ar_a = ff_grid_loglin(it_a_n, fl_a_max, fl_a_min, fl_loglin_threshold);
else
    ar_a = linspace(fl_b_bd, fl_a_max, it_a_n);
    ar_a = [0 ar_a];
    ar_a = sort(unique(ar_a));
end

%% Get Shock Grids
[ar_s, mt_z_trans, ar_stationary, ar_z] = fo_st_tauchen_jhl(fl_z_mu,fl_z_rho,fl_z_sig,it_z_n);
if (bl_display_grids)
    disp(ar_z);
    disp(mt_z_trans);
    disp(size(ar_z));
end

%% Store
armt_map = containers.Map('KeyType','char', 'ValueType','any');
armt_map('ar_a') = ar_a;
armt_map('mt_z_trans') = mt_z_trans;
armt_map('ar_stationary') = ar_stationary;
armt_map('ar_z') = ar_z;

end
