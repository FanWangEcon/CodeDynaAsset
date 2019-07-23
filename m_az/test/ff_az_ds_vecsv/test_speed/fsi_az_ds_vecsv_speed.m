%% Test Full Run Speed (Savings Distribution)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Testing the
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vecsv.html
% ff_az_ds_vecsv> program for solving the savings only dynamic
% programming problem.
%
% Here the dynamic programming problem is solved and stationary
% distribution is found. This when invoked multiple times, lead to finding
% the equilibrium prices. See
% <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3316939 Wang
% (2019)>'s Appendix for how to solve equilibrium prices in this context
% using multile cores. Specifically, if there is only one processor, the
% code could be invoked up to 12 to 15 times to find the interest rate
% using bisection. If the computer has several processors, equilibrium
% could be found using multi-section as described in Wang (2019) in 2 or 3
% iterations where the code is run concurrently across multiple processors
% within each iteration.
%
% @seealso
%
% * _SPEED_ savings only overall benchmark speed testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_speed/html/fsi_az_ds_vecsv_speed.html fsi_az_ds_vecsv_speed>
% * _PREFERENCE_ savings only preference testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_pref/html/fsi_az_ds_vecsv_pref.html fsi_az_ds_vecsv_pref>
% * _PREFERENCE_ savings only preference testing cross:
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_pref/html/fsi_az_ds_vecsv_pref_cross.html fsi_az_ds_vecsv_pref_cross>
% * _SHOCK_ savings only shock testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_shock/html/fsi_az_ds_vecsv_shock.html fsi_az_ds_vecsv_shock>
% * _SHOCK_ savings only shock testing cross:
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_shock/html/fsi_az_ds_vecsv_shock_cross.html fsi_az_ds_vecsv_shock_cross>
% * _PRICE_ savings only wage and interest rate testing cross: adjust wage and savings rate
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_price/html/fsi_az_ds_vecsv_price_cross.html fsi_az_ds_vecsv_price_cross>
% * _JOINT_ all parameters random draws joint test
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_joint/html/fsi_az_ds_vecsv_joint_rand.html fsi_az_ds_vecsv_joint_rand>
%

%% Solving the Benchmark Model

close all;
clear all;

% Start Profiling
bl_profile = true;
if (bl_profile)
    profile off;
    profile on;
end

% Start Timer
bl_time = true;
if (bl_time)
    tic;
end

% Set Parameters
bl_input_override = true;
it_param_set = 9;
[param_map, support_map] = ffs_az_set_default_param(it_param_set);
support_map('bl_time') = false;
support_map('bl_display_final_dist') = false;

% Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_get_funcgrid.html ffs_az_get_funcgrid>
[armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map, bl_input_override);

% Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html ff_az_vf_vecsv>
result_map = ff_az_vf_vecsv(param_map, support_map, armt_map, func_map);

% Call Distribution CProgram
result_map = ff_az_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

% End Timer
if (bl_time)
    toc;
end

% End Profiling
if (bl_profile)
    profile off
    profile viewer
        
    % append function name
    st_func_name = 'fsi_az_ds_vecsv_speed';
    support_map('st_profile_path') = [support_map('st_matimg_path_root') '/test/ff_az_ds_vecsv/test_speed/profile/'];    
    support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];    
    
    % support_map
    params_group = values(support_map, {'st_profile_path', ...
        'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix'});
    [st_profile_path, st_profile_prefix, st_profile_name_main, st_profile_suffix] = params_group{:};
    
    % Save
    st_file_name = [st_profile_prefix st_profile_name_main st_profile_suffix];
    profsave(profile('info'), strcat(st_profile_path, st_file_name));    
end
