%% Derive Distributions for Save + Borrow One Asset (Wrapper)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_abz_ds_wrapper(varargin)
%% FF_ABZ_DS_WRAPPER finds the stationary asset distributions
% This is a warpper function. Note that when invoked by *abz*, there are
% two shocks, one for wage, another for the borrowing interest rate.
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_set_default_param.html ffs_abz_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
% * save + borr optimized-vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
% * derive distribution f(y'(y,z)) one asset *loop*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html ff_az_ds>
% * derive distribution f(y'({y},z)) or f(y'({x,y},z)) *vectorized*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vec.html ff_az_ds_vec>
% * derive distribution f(y'({y},z)) or f(y'({x,y},z)) *semi-analytical*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vecsv.html ff_az_ds_vecsv>
%

%% Default
% # it_subset = 5 is basic invoke quick test
% # it_subset = 6 is invoke full test
% # it_subset = 7 is profiling invoke
% # it_subset = 8 is matlab publish
% # it_subset = 9 is invoke operational (only final stats) and coh graph

it_param_set = 8;
bl_input_override = true;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% Note: param_map and support_map can be adjusted here or outside to override defaults
% param_map('it_a_n') = 750;
% param_map('fl_z_r_borr_n') = 5;
% param_map('it_z_wage_n') = 15;
% param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
% param_map('fl_r_save') = 0.025;

% These parameters below for comparison with abz_fibs model
% param_map('fl_r_save') = 0.025;
% param_map('fl_r_borr') = 0.095;
param_map('fl_z_r_borr_poiss_mean') = 20;
param_map('fl_c_min') = 0.02;

% param_map('st_analytical_stationary_type') = 'loop';
% param_map('st_analytical_stationary_type') = 'vector';
param_map('st_analytical_stationary_type') = 'eigenvector';

% get armt and func map
[armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
default_params = {param_map support_map armt_map func_map};

%% Parse Parameters 1

% if varargin only has param_map and support_map,
params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};
param_map = [param_map; default_params{1}];
support_map = [support_map; default_params{2}];
if params_len >= 1 && params_len <= 2
    % If override param_map, re-generate armt and func if they are not
    % provided
    bl_input_override = true;
    [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);
else
    % Override all
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

% if profile, profile DP + Dist here
support_map('bl_profile_dist') = false;

% append function name
st_func_name = 'ff_abz_ds_wrapper';
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

bl_input_override = true;
result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);

%% Derive Distribution
% Note that the *pYisMINY* statistics show the
% proportion of households at the minimum borrowing bound, where default
% takes place, and shows therefore the fraction of households defaulting. 

if (strcmp(st_analytical_stationary_type, 'loop'))
    
    result_map = ff_az_ds(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
    
elseif (strcmp(st_analytical_stationary_type, 'vector'))
    
    result_map = ff_az_ds_vec(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
    
elseif (strcmp(st_analytical_stationary_type, 'eigenvector'))
    
    result_map = ff_az_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
    
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
