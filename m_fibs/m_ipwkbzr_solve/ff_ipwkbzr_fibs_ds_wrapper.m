%% Derive Distributions for For+Inf+Borr+Save+RShock Risky + Safe Asset Interpolated-Percentage (Wrapper)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_ipwkbzr_fibs_ds_wrapper(varargin)
%% FF_IPWKBZ_FIBS_DS_WRAPPER finds the stationary asset distributions
% This is a warpper function.

%% Default
% # it_subset = 5 is basic invoke quick test
% # it_subset = 6 is invoke full test
% # it_subset = 7 is profiling invoke
% # it_subset = 8 is matlab publish
% # it_subset = 9 is invoke operational (only final stats) and coh graph

it_param_set = 7;
[param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);

%% Change Parameter to Main Options

% Set Parameter Types
st_param_which = 'dense';

if (ismember(st_param_which, ["default"]))
    
elseif (ismember(st_param_which, ["dense"]))
    support_map('it_display_every') = 1;
    
%     param_map('it_maxiter_val') = 20;
    
    param_map('fl_coh_interp_grid_gap') = 0.05;
    param_map('it_c_interp_grid_gap') = 10^-4;
    param_map('it_w_perc_n') = 100;
    param_map('it_ak_perc_n') = param_map('it_w_perc_n');
    param_map('fl_w_interp_grid_gap') = 0.05;
    param_map('it_z_wage_n') = 5;
    param_map('it_coh_bridge_perc_n') = param_map('it_w_perc_n')/5;
    param_map('fl_z_r_infbr_n') = 11;
    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');   
    
elseif ismember(st_param_which, ["ff_ipwkbzr_ds_wrapper", "ff_ipwkbzrr_ds_wrapper"])

    param_map('fl_r_save') = 0.025;
    
    param_map('fl_r_fsv') = 0.025;
    param_map('fl_r_fbr') = 1.000;
    param_map('bl_bridge') = false;
    param_map('it_coh_bridge_perc_n') = 1;
    
    if ismember(st_param_which, ["ff_ipwkbzr_ds_wrapper"])

        % ff_ipwkbzr_evf default
        param_map('fl_z_r_borr_min') = 0.025;
        param_map('fl_z_r_borr_max') = 0.025;
        param_map('fl_z_r_borr_n') = 1;
   
    end
    
    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');    
    
end

%% Adjust Parametesr

% Note: param_map and support_map can be adjusted here or outside to override defaults
% param_map('it_w_perc_n') = 50;
% param_map('it_ak_perc_n') = param_map('it_w_perc_n');
% param_map('it_z_n') = 15;
% param_map('fl_coh_interp_grid_gap') = 0.025;
% param_map('it_c_interp_grid_gap') = 0.001;
% param_map('fl_w_interp_grid_gap') = 0.25;
% param_map('it_w_perc_n') = 100;
% param_map('it_ak_perc_n') = param_map('it_w_perc_n');
% param_map('it_z_n') = 11;
% param_map('fl_coh_interp_grid_gap') = 0.1;
% param_map('it_c_interp_grid_gap') = 10^-4;
% param_map('fl_w_interp_grid_gap') = 0.1;
% param_map('it_w_perc_n') = 100;
% param_map('fl_r_save') = 0.025;
% param_map('fl_c_min') = 0.02;

%% Set Distribution Derivation Types

% param_map('st_analytical_stationary_type') = 'loop';
% param_map('st_analytical_stationary_type') = 'vector';
param_map('st_analytical_stationary_type') = 'eigenvector';

%% Generate Grids

% get armt and func map
params_len = length(varargin);
if params_len <= 2
    [armt_map, func_map] = ffs_ipwkbzr_fibs_get_funcgrid(param_map, support_map); % 1 for override
    default_params = {param_map support_map armt_map func_map};
end

%% Parse Parameters 1

% if varargin only has param_map and support_map,
params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};
param_map = [param_map; default_params{1}];
support_map = [support_map; default_params{2}];
if params_len >= 1 && params_len <= 2
    % If override param_map, re-generate armt and func if they are not
    % provided
    [armt_map, func_map] = ffs_ipwkbzr_fibs_get_funcgrid(param_map, support_map);
else
    % Override all
    armt_map = default_params{3};
    func_map = default_params{4};
end

% if profile, profile DP + Dist here
support_map('bl_profile_dist') = false;

% append function name
st_func_name = 'ff_ipwkbzr_fibs_ds_wrapper';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters

% param_map
params_group = values(param_map, {'st_analytical_stationary_type'});
[st_analytical_stationary_type] = params_group{:};

% support_map
params_group = values(support_map, ...
    {'st_profile_path', 'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix','bl_time'});
[st_profile_path, st_profile_prefix, st_profile_name_main, st_profile_suffix, bl_time] = params_group{:};

%% Start Profiler and Timer
% Start Profile
if (it_param_set == 7)
    close all;
    profile off;
    profile on;
end

% Start Timer
if (bl_time)
    tic;
end

%% Solve DP

result_map = ff_ipwkbzr_fibs_vf_vecsv(param_map, support_map, armt_map, func_map);

%% Derive Distribution

if (strcmp(st_analytical_stationary_type, 'loop'))
    
    result_map = ff_iwkz_ds(param_map, support_map, armt_map, func_map, result_map);
    
elseif (strcmp(st_analytical_stationary_type, 'vector'))
    
    result_map = ff_iwkz_ds_vec(param_map, support_map, armt_map, func_map, result_map);
    
elseif (strcmp(st_analytical_stationary_type, 'eigenvector'))
    
    result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map);
    
end

%% End Profiler and Timer
% End Timer
if (bl_time)
    toc;
end

% End Profile
if (it_param_set == 7)
    profile off
    profile viewer
    st_file_name = [st_profile_prefix st_profile_name_main st_profile_suffix];
    profsave(profile('info'), strcat(st_profile_path, st_file_name));
end
end
