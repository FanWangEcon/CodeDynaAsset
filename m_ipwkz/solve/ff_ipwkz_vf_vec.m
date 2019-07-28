%% Risky + Safe Asset (Saving Only) Interpolated-Percentage (Vectorized)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository> 
% Table of Content.*

%%
function result_map = ff_ipwkz_vf_vec(varargin)
%% FF_IPWKZ_VF_VEC solve infinite horizon exo shock + endo asset problem
% This program solves the infinite horizon dynamic savings and risky
% capital asset problem with some shock. This is the two step solution
% with interpolation and with percentage asset grids version of
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_vf_vec.html
% ff_iwkz_vf_vec>. See that file for more descriptions. This is the
% vectorized version of the program.
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param armt_map container container with states, choices and shocks
% grids that are inputs for grid based solution algorithm
%
% @param func_map container container with function handles for
% consumption cash-on-hand etc.
%
% @return result_map container contains policy function matrix, value
% function matrix, iteration results, and policy function, value function
% and iteration results tables. 
%
% keys included in result_map:
%
% * mt_val matrix states_n by shock_n matrix of converged value function grid
% * mt_pol_a matrix states_n by shock_n matrix of converged policy function grid
% * ar_val_diff_norm array if bl_post = true it_iter_last by 1 val function
% difference between iteration
% * ar_pol_diff_norm array if bl_post = true it_iter_last by 1 policy
% function difference between iterations
% * mt_pol_perc_change matrix if bl_post = true it_iter_last by shock_n the
% proportion of grid points at which policy function changed between
% current and last iteration for each element of shock
%
% @example
%
% @include
%
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkz/paramfunc/ff_ipwkz_evf.m ff_ipwkz_evf>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkz/paramfunc/ffs_ipwkz_set_default_param.m ffs_ipwkz_set_default_param>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkz/paramfunc/ffs_ipwkz_get_funcgrid.m ffs_ipwkz_get_funcgrid>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solvepost/ff_akz_vf_post.m ff_akz_vf_post>
%

%% Default
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button

it_param_set = 4;
[param_map, support_map] = ffs_ipwkz_set_default_param(it_param_set);

% parameters can be set inside ffs_ipwkz_set_default_param or updated here
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

% get armt and func map
[armt_map, func_map] = ffs_ipwkz_get_funcgrid(param_map, support_map); % 1 for override
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
    [armt_map, func_map] = ffs_ipwkz_get_funcgrid(param_map, support_map);
else
    % Override all
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

% append function name
st_func_name = 'ff_ipwkz_vf_vec';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, {'ar_w_level', 'ar_z'});
[ar_w_level, ar_z] = params_group{:};
params_group = values(armt_map, {'ar_interp_c_grid', 'ar_interp_coh_grid', ...
    'mt_interp_coh_grid_mesh_z', 'mt_z_mesh_coh_interp_grid',...
    'mt_w_by_interp_coh_interp_grid'});
[ar_interp_c_grid, ar_interp_coh_grid, ...
    mt_interp_coh_grid_mesh_z, mt_z_mesh_coh_interp_grid, ...
    mt_w_by_interp_coh_interp_grid] = params_group{:};
params_group = values(armt_map, {'mt_coh_wkb', 'mt_z_mesh_coh_wkb'});
[mt_coh_wkb, mt_z_mesh_coh_wkb] = params_group{:};

% func_map
params_group = values(func_map, {'f_util_log', 'f_util_crra', 'f_cons'});
[f_util_log, f_util_crra, f_cons] = params_group{:};

% param_map
params_group = values(param_map, {'it_z_n', 'fl_crra', 'fl_beta', 'fl_c_min'});
[it_z_n, fl_crra, fl_beta, fl_c_min] = params_group{:};
params_group = values(param_map, {'it_maxiter_val', 'fl_tol_val', 'fl_tol_pol', 'it_tol_pol_nochange'});
[it_maxiter_val, fl_tol_val, fl_tol_pol, it_tol_pol_nochange] = params_group{:};

% support_map
params_group = values(support_map, {'bl_profile', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time', 'bl_display_defparam', 'bl_graph_evf', 'bl_display', 'it_display_every', 'bl_post'});
[bl_profile, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display_defparam, bl_graph_evf, bl_display, it_display_every, bl_post] = params_group{:};
params_group = values(support_map, {'it_display_summmat_rowmax', 'it_display_summmat_colmax'});
[it_display_summmat_rowmax, it_display_summmat_colmax] = params_group{:};

%% Initialize Output Matrixes

mt_val_cur = zeros(length(ar_interp_coh_grid),length(ar_z));
mt_val = mt_val_cur - 1;
mt_pol_a = zeros(length(ar_interp_coh_grid),length(ar_z));
mt_pol_a_cur = mt_pol_a - 1;
mt_pol_k = zeros(length(ar_interp_coh_grid),length(ar_z));
mt_pol_k_cur = mt_pol_k - 1;
mt_pol_idx = zeros(length(ar_interp_coh_grid),length(ar_z));

%% Initialize Convergence Conditions

bl_vfi_continue = true;
it_iter = 0;
ar_val_diff_norm = zeros([it_maxiter_val, 1]);
ar_pol_diff_norm = zeros([it_maxiter_val, 1]);
mt_pol_perc_change = zeros([it_maxiter_val, it_z_n]);

%% Pre-calculate u(c)
% Interpolation, see
% <https://fanwangecon.github.io/M4Econ/support/speed/partupdate/fs_u_c_partrepeat_main.html
% fs_u_c_partrepeat_main> for why interpolate over u(c)

% Evaluate
if (fl_crra == 1)
    ar_interp_u_of_c_grid = f_util_log(ar_interp_c_grid);
    fl_u_neg_c = f_util_log(fl_c_min);
else
    ar_interp_u_of_c_grid = f_util_crra(ar_interp_c_grid);
    fl_u_neg_c = f_util_crra(fl_c_min);
end
ar_interp_u_of_c_grid(ar_interp_c_grid <= fl_c_min) = fl_u_neg_c;

% Get Interpolant
f_grid_interpolant_spln = griddedInterpolant(ar_interp_c_grid, ar_interp_u_of_c_grid, 'spline', 'nearest');

%% Iterate Value Function
% Loop solution with 4 nested loops
%
% # loop 1: over exogenous states
% # loop 2: over endogenous states
% # loop 3: over choices
% # loop 4: add future utility, integration--loop over future shocks
%

% Start Profile
if (bl_profile)
    close all;
    profile off;
    profile on;
end

% Start Timer
if (bl_time)
    tic;
end

% Value Function Iteration
while bl_vfi_continue
    it_iter = it_iter + 1;
    
    
    %% Interpolate (1) reacahble v(coh(k(w,z),b(w,z),z),z) given v(coh, z)
    % v(coh,z) solved on ar_interp_coh_grid, ar_z grids, see
    % ffs_ipwkz_get_funcgrid.m. Generate interpolant based on that, Then
    % interpolate for the coh reachable levels given the k(w,z) percentage
    % choice grids in the second stage of the problem

    % Generate Interpolant for v(coh,z)
    f_grid_interpolant_value = griddedInterpolant(...
        mt_z_mesh_coh_interp_grid', mt_interp_coh_grid_mesh_z', mt_val_cur', 'linear', 'nearest');
    
    % Interpoalte for v(coh(k(w,z),b(w,z),z),z)
    mt_val_wkb_interpolated = f_grid_interpolant_value(mt_z_mesh_coh_wkb, mt_coh_wkb);
    
    %% Solve Second Stage Problem k*(w,z)
    % This is the key difference between this function and
    % <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_functions.html
    % ffs_akz_set_functions> which solves the two stages jointly    
    % Interpolation first, because solution coh grid is not the same as all
    % points reachable by k and b choices given w. 
    
    support_map('bl_graph_evf') = false;
    if (it_iter == (it_maxiter_val + 1))
        support_map('bl_graph_evf') = bl_graph_evf;
    end
    
    [mt_ev_condi_z_max, ~, mt_ev_condi_z_max_kp, ~] = ...
        ff_ipwkz_evf(mt_val_wkb_interpolated, param_map, support_map, armt_map);       
    
    %% Solve First Stage Problem w*(z) given k*(w,z)
       
    % loop 1: over exogenous states
    for it_z_i = 1:length(ar_z)

        % State Array fixed        
        ar_coh_z = mt_interp_coh_grid_mesh_z(:,it_z_i);

        % Get 2nd Stage Arrays
        ar_ev_condi_z_max_z = mt_ev_condi_z_max(:, it_z_i);        
        ar_w_level_kstar_z = mt_ev_condi_z_max_kp(:, it_z_i);

        % Interpolate (2) k*(ar_w_perc) from k*(ar_w_level,z)
        % There are two w=k'+b' arrays. ar_w_level is the level even grid based
        % on which we solve the 2nd stage problem in ff_ipwkz_evf.m. Here for
        % each coh level, we have a different vector of w levels, but the same
        % vector of percentage ws. So we need to interpolate to get the optimal
        % k* and b* choices at each percentage level of w. 
        f_interpolante_w_level_kstar_z = griddedInterpolant(ar_w_level, ar_w_level_kstar_z', 'linear', 'nearest');
        
        % Interpolant for (3) EV(k*(ar_w_perc),Z)
        f_interpolante_ev_condi_z_max_z = griddedInterpolant(ar_w_level, ar_ev_condi_z_max_z', 'linear', 'nearest');
        
        % Interpolat (2) and (3) 
        mt_w_kstar_interp_z = f_interpolante_w_level_kstar_z(mt_w_by_interp_coh_interp_grid);
        mt_w_astar_interp_z = mt_w_by_interp_coh_interp_grid - mt_w_kstar_interp_z;
        mt_ev_condi_z_max_interp_z = f_interpolante_ev_condi_z_max_z(mt_w_by_interp_coh_interp_grid);
        
        % Consumption
        % Note that compared to
        % <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_functions.html
        % ffs_akz_set_functions> the mt_c here is much smaller the same
        % number of columns (states) as in the ffs_akz_set_functions file,
        % but the number of rows equal to ar_w length.
        mt_c = f_cons(ar_coh_z', mt_w_astar_interp_z, mt_w_kstar_interp_z);
        
        % Interpolate for (4) EVAL current utility: N by N, f_util defined earlier
        mt_utility = f_grid_interpolant_spln(mt_c);
                
        % EVAL add on future utility, N by N + N by N
        % previously mt_ev_condi_z_max_interp_z was N by 1, not N by N
        mt_utility = mt_utility + fl_beta*mt_ev_condi_z_max_interp_z;

        % Optimization: remember matlab is column major, rows must be
        % choices, columns must be states
        % <https://en.wikipedia.org/wiki/Row-_and_column-major_order COLUMN-MAJOR>
        [ar_opti_val1_z, ar_opti_idx_z] = max(mt_utility);
        mt_val(:,it_z_i) = ar_opti_val1_z;
        
        % Generate Linear Opti Index
        [it_choies_n, it_states_n] = size(mt_utility);
        ar_add_grid = linspace(0, it_choies_n*(it_states_n-1), it_states_n);
        ar_opti_linear_idx_z = ar_opti_idx_z + ar_add_grid;

        % Store Optimal
        mt_pol_a(:,it_z_i) = mt_w_astar_interp_z(ar_opti_linear_idx_z);
        mt_pol_k(:,it_z_i) = mt_w_kstar_interp_z(ar_opti_linear_idx_z);        
        if (it_iter == (it_maxiter_val + 1))
            mt_pol_idx(:,it_z_i) = ar_opti_linear_idx_z;
        end

    end
    
    %% Check Tolerance and Continuation
    
    % Difference across iterations
    ar_val_diff_norm(it_iter) = norm(mt_val - mt_val_cur);
    ar_pol_diff_norm(it_iter) = norm(mt_pol_a - mt_pol_a_cur) + norm(mt_pol_k - mt_pol_k_cur);
    ar_pol_a_perc_change = sum((mt_pol_a ~= mt_pol_a_cur))/(length(ar_interp_coh_grid));
    ar_pol_k_perc_change = sum((mt_pol_k ~= mt_pol_k_cur))/(length(ar_interp_coh_grid));    
    mt_pol_perc_change(it_iter, :) = mean([ar_pol_a_perc_change;ar_pol_k_perc_change]);
    
    % Update
    mt_val_cur = mt_val;
    mt_pol_a_cur = mt_pol_a;
    mt_pol_k_cur = mt_pol_k;
    
    % Print Iteration Results
    if (bl_display && (rem(it_iter, it_display_every)==0))
        fprintf('VAL it_iter:%d, fl_diff:%d, fl_diff_pol:%d\n', ...
            it_iter, ar_val_diff_norm(it_iter), ar_pol_diff_norm(it_iter));
        tb_valpol_iter = array2table([mean(mt_val_cur,1);...
                                      mean(mt_pol_a_cur,1); ...
                                      mean(mt_pol_k_cur,1); ...
                                      mt_val_cur(length(ar_interp_coh_grid),:); ...
                                      mt_pol_a_cur(length(ar_interp_coh_grid),:); ...
                                      mt_pol_k_cur(length(ar_interp_coh_grid),:)]);
        tb_valpol_iter.Properties.VariableNames = strcat('z', string((1:size(mt_val_cur,2))));
        tb_valpol_iter.Properties.RowNames = {'mval', 'map', 'mak', 'Hval', 'Hap', 'Hak'};
        disp('mval = mean(mt_val_cur,1), average value over a')
        disp('map  = mean(mt_pol_a_cur,1), average choice over a')
        disp('mkp  = mean(mt_pol_k_cur,1), average choice over k')
        disp('Hval = mt_val_cur(it_ameshk_n,:), highest a state val')
        disp('Hap = mt_pol_a_cur(it_ameshk_n,:), highest a state choice')
        disp('mak = mt_pol_k_cur(it_ameshk_n,:), highest k state choice')                
        disp(tb_valpol_iter);
    end
    
    % Continuation Conditions:
    % 1. if value function convergence criteria reached
    % 2. if policy function variation over iterations is less than
    % threshold
    if (it_iter == (it_maxiter_val + 1))
        bl_vfi_continue = false;
    elseif ((it_iter == it_maxiter_val) || ...
            (ar_val_diff_norm(it_iter) < fl_tol_val) || ...
            (sum(ar_pol_diff_norm(max(1, it_iter-it_tol_pol_nochange):it_iter)) < fl_tol_pol))
        % Fix to max, run again to save results if needed
        it_iter_last = it_iter;
        it_iter = it_maxiter_val;        
    end
    
end

% End Timer
if (bl_time)
    toc;
end

% End Profile
if (bl_profile)
    profile off
    profile viewer
    st_file_name = [st_profile_prefix st_profile_name_main st_profile_suffix];
    profsave(profile('info'), strcat(st_profile_path, st_file_name));
end

%% Process Optimal Choices

result_map = containers.Map('KeyType','char', 'ValueType','any');
result_map('mt_val') = mt_val;
result_map('mt_pol_idx') = mt_pol_idx;

result_map('cl_mt_pol_coh') = {mt_interp_coh_grid_mesh_z, zeros(1)};
result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
result_map('cl_mt_pol_k') = {mt_pol_k, zeros(1)};
result_map('cl_mt_pol_c') = {f_cons(mt_interp_coh_grid_mesh_z, mt_pol_a, mt_pol_k), zeros(1)};
result_map('ar_st_pol_names') = ["cl_mt_pol_coh", "cl_mt_pol_a", "cl_mt_pol_k", "cl_mt_pol_c"];

if (bl_post)
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);

    % graphing based on coh_wkb, but that does not match optimal choice
    % matrixes for graphs.
    armt_map('mt_coh_wkb') = mt_interp_coh_grid_mesh_z;
    armt_map('it_ameshk_n') = length(ar_interp_coh_grid);
    armt_map('ar_a_meshk') = mt_interp_coh_grid_mesh_z(:,1);
    armt_map('ar_k_mesha') = zeros(size(mt_interp_coh_grid_mesh_z(:,1)) + 0);

    result_map = ff_akz_vf_post(param_map, support_map, armt_map, func_map, result_map);
end

%% Display Various Containers

if (bl_display_defparam)

    %% Display 1 support_map
    fft_container_map_display(support_map, it_display_summmat_rowmax, it_display_summmat_colmax);

    %% Display 2 armt_map
    fft_container_map_display(armt_map, it_display_summmat_rowmax, it_display_summmat_colmax);

    %% Display 3 param_map
    fft_container_map_display(param_map, it_display_summmat_rowmax, it_display_summmat_colmax);

    %% Display 4 func_map
    fft_container_map_display(func_map, it_display_summmat_rowmax, it_display_summmat_colmax);

    %% Display 5 result_map
    fft_container_map_display(result_map, it_display_summmat_rowmax, it_display_summmat_colmax);

end

end
