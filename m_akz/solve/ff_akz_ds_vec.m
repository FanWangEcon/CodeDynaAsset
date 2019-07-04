%% Derive Two Asset (Risky + Safe) and Choices/Outcomes Distribution (Vectorized)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% This uses
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vec.html
% ff_az_ds_vec>, which works for single and multiple assets.
%
% The function here works with both
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vecsv.html
% ff_akz_vf_vecsv> as well as
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_wkz_vf_vecsv.html
% ff_wkz_vf_vecsv>. Results are identical, but _ff_wkz_vf_vecsv_ is
% significantly faster.
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vecsv.html ff_az_vf_vecsv>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_wkz_vf_vecsv.html ff_wkz_vf_vecsv>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vec.html ff_az_ds_vec>
%
% @seealso
%
% * derive distribution f(y'(y,z)) one asset *loop*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html ff_az_ds>
% * derive distribution f(y'({x,y},z)) two assets *loop*: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_ds.html ff_akz_ds>
% * derive distribution f(y'({x,y},z, *z'*)) two assets *loop*: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_ds.html ff_iwkz_ds>
% * derive distribution f(y'({y},z)) or f(y'({x,y},z)) *vectorized*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vec.html ff_az_ds_vec>
% * derive distribution f(y'({y},z, *z'*)) or f(y'({x,y},z, *z'*)) *vectorized*: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_ds_vec.html ff_iwkz_ds_vec>
% * derive distribution f(y'({y},z)) or f(y'({x,y},z)) *semi-analytical*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vecsv.html ff_az_ds_vecsv>
% * derive distribution f(y'({y},z, *z'*)) or f(y'({x,y},z, *z'*)) *semi-analytical*: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_ds_vecsv.html ff_iwkz_ds_vecsv>
%

%% Set Parameters Main
% Options for Distribution solutions
%
% # it_subset = 5 is basic invoke quick test
% # it_subset = 6 is invoke full test
% # it_subset = 7 is profiling invoke
% # it_subset = 8 is matlab publish
% # it_subset = 9 is invoke operational (only final stats) and coh graph
%

close all;
clear all;

it_param_set = 8;
st_akz_or_wkz = 'wkz';

%% Get Parameters
% Note that akz and wkz share the smae funcgrid and default_param functions

% Set Parameters
bl_input_override = true;
[param_map, support_map] = ffs_akz_set_default_param(it_param_set);
support_map('bl_profile_dist') = false;

% Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html ffs_akz_get_funcgrid>
[armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override

%% Alternative 1: Solving the Dynamic Programming Problem (AKZ)

if (strcmp(st_akz_or_wkz, 'akz'))
    result_map = ff_akz_vf_vecsv(param_map, support_map, armt_map, func_map);
end

%% Alternative 2: Solving the Dynamic Programming Problem (AWZ)

if (strcmp(st_akz_or_wkz, 'wkz'))
    result_map = ff_wkz_vf_vecsv(param_map, support_map, armt_map, func_map);
end

%% Distribution Start Profiler and Timer

% Start Profiling
if (it_param_set == 7)
    profile off;
    profile on;
end

% Start Timer
if (support_map('bl_time'))
    tic;
end

% Get Profile name
st_profile_name_main = support_map('st_profile_name_main');

%% Derive Distribution

% Call Distribution Program <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vec.html ff_az_ds_vec>
result_map = ff_az_ds_vec(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

%% End profiler and Timer

% End Timer
if (support_map('bl_time'))
    toc;
end

% End Profiling
if (it_param_set == 7)
    profile off
    profile viewer
        
    % append function name
    st_func_name = 'ff_akz_ds_vec';
    support_map('st_profile_path') = [support_map('st_matimg_path_root') '/solve/profile/'];    
    support_map('st_profile_name_main') = [st_func_name st_profile_name_main];    
    
    % support_map
    params_group = values(support_map, {'st_profile_path', ...
        'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix'});
    [st_profile_path, st_profile_prefix, st_profile_name_main, st_profile_suffix] = params_group{:};
    
    % Save
    st_file_name = [st_profile_prefix st_profile_name_main st_profile_suffix];
    profsave(profile('info'), strcat(st_profile_path, st_file_name));    
end
