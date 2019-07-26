%% Test Full Run Speed (Risky + Safe Asets + Interpolated)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_test_analyze.html ff_az_test_analyze>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_akz_set_default_param.html ffs_akz_set_default_param>
%
% @seealso
%
% * _SPEED_ risky + safe (two-step interp) overall benchmark speed testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_speed/html/fsi_ikwz_ds_vecsv_speed.html fsi_ikwz_ds_vecsv_speed>
% * _PREFERENCE_ risky + safe (two-step interp) preference testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_pref/html/fsi_ikwz_ds_vecsv_pref.html fsi_ikwz_ds_vecsv_pref>
% * _PREFERENCE_ risky + safe (two-step interp) preference testing *cross*:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_pref/html/fsi_ikwz_ds_vecsv_pref_cross.html fsi_ikwz_ds_vecsv_pref_cross>
% * _SHOCK_ risky + safe (two-step interp) shock testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_shock/html/fsi_ikwz_ds_vecsv_shock.html fsi_ikwz_ds_vecsv_shock>
% * _SHOCK_ risky + safe (two-step interp) shock testing *cross*:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_shock/html/fsi_ikwz_ds_vecsv_shock_cross.html fsi_ikwz_ds_vecsv_shock_cross>
% * _CAPITAL_ risky + safe (two-step interp) capital return testing: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_prod/html/fsi_ikwz_ds_vecsv_prod.html fsi_ikwz_ds_vecsv_prod>
% * _CAPITAL_ risky + safe (two-step interp) capital return testing *cross*:  <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_prod/html/fsi_ikwz_ds_vecsv_prod_cross.html fsi_ikwz_ds_vecsv_prod_cross>
% * _PRICE_ risky + safe (two-step interp) wage and interest rate testing *cross*: adjust wage and savings rate
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_price/html/fsi_ikwz_ds_vecsv_price_cross.html fsi_ikwz_ds_vecsv_price_cross>
% * _JOINT_ all parameters random draws joint test
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/test/ff_ikwz_ds_vecsv/test_joint/html/fsi_ikwz_ds_vecsv_joint_rand.html fsi_ikwz_ds_vecsv_joint_rand>
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
it_param_set = 9;
[param_map, support_map] = ffs_akz_set_default_param(it_param_set);
support_map('bl_time') = false;
support_map('bl_display_final_dist') = false;

% Call Grid Generator <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_akz_get_funcgrid.html ffs_akz_get_funcgrid>
[armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map);

% Call Dynamic Programming Problem <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_iwkz_vf_vecsv.html ff_iwkz_vf_vecsv>
result_map = ff_iwkz_vf_vecsv(param_map, support_map, armt_map, func_map);

% Call Distribution CProgram
result_map = ff_iwkz_ds_vec(param_map, support_map, armt_map, func_map, result_map);

% End Timer
if (bl_time)
    toc;
end

% End Profiling
if (bl_profile)
    profile off
    profile viewer

    % append function name
    st_func_name = 'fsi_iwkz_ds_vecsv_speed';
    support_map('st_profile_path') = [support_map('st_matimg_path_root') '/test/ff_ikwz_ds_vecsv/test_speed/profile/'];
    support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];

    % support_map
    params_group = values(support_map, {'st_profile_path', ...
        'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix'});
    [st_profile_path, st_profile_prefix, st_profile_name_main, st_profile_suffix] = params_group{:};

    % Save
    st_file_name = [st_profile_prefix st_profile_name_main st_profile_suffix];
    profsave(profile('info'), strcat(st_profile_path, st_file_name));
end
