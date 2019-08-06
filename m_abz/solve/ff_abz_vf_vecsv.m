%% Solve Save + Borr (RShock) Dynamic Programming Problem (Optimized-Vectorized)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function result_map = ff_abz_vf_vecsv(varargin)
%% FF_ABZ_VF_VECSV solve infinite horizon exo shock + endo asset problem
% This program solves the infinite horizon dynamic single asset and two
% shocks problem with vectorized codes.
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf.html
% ff_abz_vf> shows looped codes.
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf_vec.html
% ff_abz_vf_vec> shows vectorized codes. This file shows vectorized codes
% that is faster but is more memory intensive.
%
% The borrowing problem is similar to the savings problem. The main
% addition here in comparison to the savings only code
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html
% ff_az_vf_vec> is the ability to deal with default, as well as an
% additional shock to the borrowing interest rate.
%
% See
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf_vec.html
% ff_az_vf_vec> how vectorization works within this structure.
%
% This _optimized-vectorized_ solution method provides very large speed
% improvements for this infinite horizon problem because the u(c(z,a,a'))
% calculation within each iteration is identical. Generally the idea is to
% identify inside iteration whether the model is infinite horizon or
% life-cycle based where repeat calculations are taking place. If such
% calculations can be identified, then potentially they could be stored and
% retrieved during future iterations/periods rather than recomputed every
% time. This saves time.
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
%    % Get Default Parameters
%    it_param_set = 2;
%    [param_map, support_map] = ffs_abz_set_default_param(it_param_set);
%    % Chnage param_map keys for borrowing
%    param_map('fl_b_bd') = -20; % borrow bound
%    param_map('bl_default') = false; % true if allow for default
%    param_map('fl_c_min') = 0.0001; % u(c_min) when default
%    % Change Keys in param_map
%    param_map('it_a_n') = 500;
%    param_map('fl_z_r_borr_n') = 5;
%    param_map('it_z_wage_n') = 15;
%    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
%    param_map('fl_a_max') = 100;
%    param_map('fl_w') = 1.3;
%    % Change Keys support_map
%    support_map('bl_display') = false;
%    support_map('bl_post') = true;
%    support_map('bl_display_final') = false;
%    % Call Program with external parameters that override defaults.
%    ff_abz_vf_vecsv(param_map, support_map);
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_set_default_param.html ffs_abz_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_vf_post.html ff_az_vf_post>
%
% @seealso
%
% * save loop: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf.html ff_az_vf>
% * save vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vec.html ff_az_vf_vec>
% * save optimized-vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html ff_az_vf_vecsv>
% * save + borr loop: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf.html ff_abz_vf>
% * save + borr vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf_vec.html ff_abz_vf_vec>
% * save + borr optimized-vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
%

%% Default
%
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button
%

it_param_set = 4;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% Note: param_map and support_map can be adjusted here or outside to override defaults
% param_map('it_a_n') = 750;
% param_map('fl_z_r_borr_n') = 5;
% param_map('it_z_wage_n') = 15;
% param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
% param_map('fl_r_save') = 0.025;
% param_map('fl_z_r_borr_poiss_mean') = 1.75;

[armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map); % 1 for override
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
    [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map);
else
    % Override all
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

% append function name
st_func_name = 'ff_abz_vf_vecsv';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, {'ar_a', 'mt_z_trans', 'ar_z_r_borr_mesh_wage', 'ar_z_wage_mesh_r_borr'});
[ar_a, mt_z_trans, ar_z_r_borr_mesh_wage, ar_z_wage_mesh_r_borr] = params_group{:};

% func_map
params_group = values(func_map, {'f_util_log', 'f_util_crra', 'f_cons_checkcmin', 'f_coh', 'f_cons_coh'});
[f_util_log, f_util_crra, f_cons_checkcmin, f_coh, f_cons_coh] = params_group{:};

% param_map
params_group = values(param_map, {'it_a_n', 'it_z_n', 'fl_crra', 'fl_beta', 'fl_c_min',...
    'fl_nan_replace', 'bl_default', 'fl_default_aprime'});
[it_a_n, it_z_n, fl_crra, fl_beta, fl_c_min, ...
    fl_nan_replace, bl_default, fl_default_aprime] = params_group{:};
params_group = values(param_map, {'it_maxiter_val', 'fl_tol_val', 'fl_tol_pol', 'it_tol_pol_nochange'});
[it_maxiter_val, fl_tol_val, fl_tol_pol, it_tol_pol_nochange] = params_group{:};

% support_map
params_group = values(support_map, {'bl_profile', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time', 'bl_display_defparam', 'bl_display', 'it_display_every', 'bl_post'});
[bl_profile, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display_defparam, bl_display, it_display_every, bl_post] = params_group{:};
params_group = values(support_map, {'it_display_summmat_rowmax', 'it_display_summmat_colmax'});
[it_display_summmat_rowmax, it_display_summmat_colmax] = params_group{:};

%% Initialize Output Matrixes
% include mt_pol_idx which we did not have in looped code

mt_val_cur = zeros(it_a_n, it_z_n);
mt_val = mt_val_cur - 1;
mt_pol_a = zeros(it_a_n, it_z_n);
mt_pol_a_cur = mt_pol_a - 1;
mt_pol_idx = zeros(it_a_n, it_z_n);

% We did not need these in ff_abz_vf or ff_abz_vf_vec
% see
% <https://fanwangecon.github.io/M4Econ/support/speed/partupdate/fs_u_c_partrepeat_main.html
% fs_u_c_partrepeat_main> for why store using cells.
cl_u_c_store = cell([it_z_n, 1]);
cl_c_valid_idx = cell([it_z_n, 1]);

%% Initialize Convergence Conditions

bl_vfi_continue = true;
it_iter = 0;
ar_val_diff_norm = zeros([it_maxiter_val, 1]);
ar_pol_diff_norm = zeros([it_maxiter_val, 1]);
mt_pol_perc_change = zeros([it_maxiter_val, it_z_n]);

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

    %% Solve Optimization Problem Current Iteration
    % Only this segment of code differs between ff_abz_vf and ff_abz_vf_vec
    % Store in cells results and retrieve, this is more memory intensive
    % than ff_abz_vf_vec.

    % loop 1: over exogenous states
    for it_z_i = 1:it_z_n

        % Current Shock
        fl_z_r_borr = ar_z_r_borr_mesh_wage(it_z_i);
        fl_z_wage = ar_z_wage_mesh_r_borr(it_z_i);

        % cash-on-hand
        ar_coh = f_coh(fl_z_r_borr, fl_z_wage, ar_a);

        % Consumption and u(c) only need to be evaluated once
        if (it_iter == 1)

            % Consumption: fl_z = 1 by 1, ar_a = 1 by N, ar_a' = N by 1
            % mt_c is N by N: matrix broadcasting, expand to matrix from arrays
            mt_c = f_cons_coh(ar_coh, ar_a');

            % EVAL current utility: N by N, f_util defined earlier
            % slightly faster to explicitly write function
            if (fl_crra == 1)
                mt_utility = log(mt_c);
                fl_u_cmin = f_util_log(fl_c_min);
            else
                % slightly faster if write function here directly, but
                % speed gain is very small, more important to have single
                % location control of functions.
                mt_utility = f_util_crra(mt_c);
                fl_u_cmin = f_util_crra(fl_c_min);
            end

            % Eliminate Complex Numbers
            mt_it_c_valid_idx = (mt_c <= fl_c_min);
            mt_utility(mt_it_c_valid_idx) = fl_u_cmin;

            % Store in cells
            cl_u_c_store{it_z_i} = mt_utility;
            cl_c_valid_idx{it_z_i} = mt_it_c_valid_idx;

        end

        % f(z'|z)
        ar_z_trans_condi = mt_z_trans(it_z_i,:);

        % EVAL EV((A',K'),Z'|Z) = V((A',K'),Z') x p(z'|z)', (N by Z) x (Z by 1) = N by 1
        % Note: transpose ar_z_trans_condi from 1 by Z to Z by 1
        % Note: matrix multiply not dot multiply
        mt_evzp_condi_z = mt_val_cur * ar_z_trans_condi';

        % EVAL add on future utility, N by N + N by 1
        mt_utility = cl_u_c_store{it_z_i} + fl_beta*mt_evzp_condi_z;

        % Index update
        % using the method below is much faster than index replace
        % see <https://fanwangecon.github.io/M4Econ/support/speed/index/fs_subscript.html fs_subscript>
        mt_it_c_valid_idx = cl_c_valid_idx{it_z_i};
        % Default or Not Utility Handling
        if (bl_default)
            % if default: only today u(cmin), transition out next period, debt wiped out
            fl_v_default = fl_u_cmin + fl_beta*mt_evzp_condi_z(ar_a == fl_default_aprime);
            mt_utility = mt_utility.*(~mt_it_c_valid_idx) + fl_v_default*(mt_it_c_valid_idx);
        else
            % if default is not allowed: v = u(cmin)
            mt_utility = mt_utility.*(~mt_it_c_valid_idx) + fl_nan_replace*(mt_it_c_valid_idx);
        end

        % Optimization: remember matlab is column major, rows must be
        % choices, columns must be states
        % <https://en.wikipedia.org/wiki/Row-_and_column-major_order COLUMN-MAJOR>
        % mt_utility is N by N, rows are choices, cols are states.
        [ar_opti_val_z, ar_opti_idx_z] = max(mt_utility);
        ar_opti_aprime_z = ar_a(ar_opti_idx_z);
        ar_opti_c_z = f_cons_coh(ar_coh, ar_opti_aprime_z);

        % Handle Default is optimal or not
        if (bl_default)
            % if defaulting is optimal choice, at these states, not required
            % to default, non-default possible, but default could be optimal
            ar_opti_aprime_z(ar_opti_c_z <= fl_c_min) = fl_default_aprime;
            ar_opti_idx_z(ar_opti_c_z <= fl_c_min) = find(ar_a == fl_default_aprime);
        else
            % if default is not allowed, then next period same state as now
            % this is absorbing state, this is the limiting case, single
            % state space point, lowest a and lowest shock has this.
            ar_opti_aprime_z(ar_opti_c_z <= fl_c_min) = ar_a(ar_opti_c_z <= fl_c_min);
        end

        % store optimal values
        mt_val(:,it_z_i) = ar_opti_val_z;
        mt_pol_a(:,it_z_i) = ar_opti_aprime_z;

        if (it_iter == (it_maxiter_val + 1))
            mt_pol_idx(:,it_z_i) = ar_opti_idx_z;
        end
    end

    %% Check Tolerance and Continuation

    % Difference across iterations
    ar_val_diff_norm(it_iter) = norm(mt_val - mt_val_cur);
    ar_pol_diff_norm(it_iter) = norm(mt_pol_a - mt_pol_a_cur);
    mt_pol_perc_change(it_iter, :) = sum((mt_pol_a ~= mt_pol_a_cur))/(it_a_n);

    % Update
    mt_val_cur = mt_val;
    mt_pol_a_cur = mt_pol_a;

    % Print Iteration Results
    if (bl_display && (rem(it_iter, it_display_every)==0))
        fprintf('VAL it_iter:%d, fl_diff:%d, fl_diff_pol:%d\n', ...
            it_iter, ar_val_diff_norm(it_iter), ar_pol_diff_norm(it_iter));
        tb_valpol_iter = array2table([mean(mt_val_cur,1); mean(mt_pol_a_cur,1); ...
            mt_val_cur(it_a_n,:); mt_pol_a_cur(it_a_n,:)]);
        tb_valpol_iter.Properties.VariableNames = strcat('z', string((1:size(mt_val_cur,2))));
        tb_valpol_iter.Properties.RowNames = {'mval', 'map', 'Hval', 'Hap'};
        disp('mval = mean(mt_val_cur,1), average value over a')
        disp('map  = mean(mt_pol_a_cur,1), average choice over a')
        disp('Hval = mt_val_cur(it_a_n,:), highest a state val')
        disp('Hap = mt_pol_a_cur(it_a_n,:), highest a state choice')
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

result_map('cl_mt_val') = {mt_val, zeros(1)};
result_map('cl_mt_coh') = {f_coh(ar_z_r_borr_mesh_wage, ar_z_wage_mesh_r_borr, ar_a'), zeros(1)};
result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
result_map('cl_mt_pol_c') = {f_cons_checkcmin(ar_z_r_borr_mesh_wage, ar_z_wage_mesh_r_borr, ar_a', mt_pol_a), zeros(1)};
result_map('ar_st_pol_names') = ["cl_mt_val", "cl_mt_pol_a", "cl_mt_coh", "cl_mt_pol_c"];

if (bl_post)
    bl_input_override = true;
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);
    result_map = ff_az_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
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