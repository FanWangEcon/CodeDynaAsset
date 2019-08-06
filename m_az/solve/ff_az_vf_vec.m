%% Solve One Asset Dynamic Programming Problem (Vectorized)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function result_map = ff_az_vf_vec(varargin)
%% FF_AZ_VF_VEC solve infinite horizon exo shock + endo asset problem
% This program solves the infinite horizon dynamic single asset and single
% shock problem with vectorized codes.
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf.html
% ff_az_vf> shows looped codes. The solution is the same.
%
% The vectorization takes advantage of implicit parallization that modern
% computers have when <https://en.wikipedia.org/wiki/SIMD same instructions
% are given for different blocks of data> With vectorization, we face a
% tradeoff between memory and speed. Suppose we have many shock points and
% many states points, if we build all states and choices into one single
% matrix and compute consumption, utility, etc over that entire matrix, that
% might be more efficient than computing consumption, utility, etc by
% subset of that matrix over a loop, but there is time required for
% generating that large input matrix, and if there are too many states, a
% computer could run out of memory.
%
% The design philosophy here is that we vectorize the endogenous states and
% choices into matrixes, but do not include the exogeous states (shocks).
% The exogenous shocks remain looped. This means we can potentially have
% multiple shock variables discretized over a large number of shock states,
% and the computer would not run into memory problems. The speed gain from
% vectoring the rest of the problem conditional on shocks is very large
% compared to the pure looped version of the problem. Even if more memory
% is available, including the exogenous states in the vectorization process
% might not be speed improving.
%
% Note one key issue is whether a programming language is
% <https://en.wikipedia.org/wiki/Row-_and_column-major_order row or column
% major> depending on which, states should be rows or columns.
%
% Another programming issue is the idea of *broadcasting* vs matrix
% algebra, both are used here. Since Matlab R2016b,
% <https://blogs.mathworks.com/loren/2016/10/24/matlab-arithmetic-expands-in-r2016b/
% matrix broadcasting> has been allowed, which means the sum of a N by 1
% and 1 by M is N by M. This is unrelated to matrix algebra. Matrix array
% broadcasting is very useful because it reduces the dimensionality of our
% model input state and choice and shock vectors, offering greater code
% clarity.
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
%    % Change Keys in param_map
%    param_map('it_a_n') = 500;
%    param_map('it_z_n') = 11;
%    param_map('fl_a_max') = 100;
%    param_map('fl_w') = 1.3;
%    % Change Keys support_map
%    support_map('bl_display') = false;
%    support_map('bl_post') = true;
%    support_map('bl_display_final') = false;
%    % Call Program with external parameters that override defaults.
%    ff_az_vf_vec(param_map, support_map);
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html ffs_az_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_get_funcgrid.html ffs_az_get_funcgrid>
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
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button

it_param_set = 4;
bl_input_override = true;
[param_map, support_map] = ffs_az_set_default_param(it_param_set);

% Note: param_map and support_map can be adjusted here or outside to override defaults
% param_map('it_a_n') = 750;
% param_map('it_z_n') = 15;

[armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
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
    [armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map, bl_input_override);
else
    % Override all
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

% append function name
st_func_name = 'ff_az_vf_vec';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, {'ar_a', 'mt_z_trans', 'ar_z'});
[ar_a, mt_z_trans, ar_z] = params_group{:};
% func_map
params_group = values(func_map, {'f_util_log', 'f_util_crra', 'f_cons', 'f_coh'});
[f_util_log, f_util_crra, f_cons, f_coh] = params_group{:};
% param_map
params_group = values(param_map, {'it_a_n', 'it_z_n', 'fl_crra', 'fl_beta', 'fl_nan_replace'});
[it_a_n, it_z_n, fl_crra, fl_beta, fl_nan_replace] = params_group{:};
params_group = values(param_map, {'it_maxiter_val', 'fl_tol_val', 'fl_tol_pol', 'it_tol_pol_nochange'});
[it_maxiter_val, fl_tol_val, fl_tol_pol, it_tol_pol_nochange] = params_group{:};

% support_map
params_group = values(support_map, {'bl_profile', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time', 'bl_display', 'it_display_every', 'bl_post'});
[bl_profile, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display, it_display_every, bl_post] = params_group{:};

%% Initialize Output Matrixes
% include mt_pol_idx which we did not have in looped code

mt_val_cur = zeros(length(ar_a),length(ar_z));
mt_val = mt_val_cur - 1;
mt_pol_a = zeros(length(ar_a),length(ar_z));
mt_pol_a_cur = mt_pol_a - 1;
mt_pol_idx = zeros(length(ar_a),length(ar_z));

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
    % Only this segment of code differs between ff_az_vf and ff_az_vf_vec

    % loop 1: over exogenous states
    % incorporating these shocks into vectorization has high memory burden
    % but insignificant speed gains. Keeping this loop allows for large
    % number of shocks without overwhelming memory
    for it_z_i = 1:length(ar_z)

        % Current Shock
        fl_z = ar_z(it_z_i);

        % Consumption: fl_z = 1 by 1, ar_a = 1 by N, ar_a' = N by 1
        % mt_c is N by N: matrix broadcasting, expand to matrix from arrays
        mt_c = f_cons(fl_z, ar_a, ar_a');

        % EVAL current utility: N by N, f_util defined earlier
        if (fl_crra == 1)
            mt_utility = f_util_log(mt_c);
        else
            mt_utility = f_util_crra(mt_c);
        end

        % f(z'|z), 1 by Z
        ar_z_trans_condi = mt_z_trans(it_z_i,:);

        % EVAL EV((A',K'),Z'|Z) = V((A',K'),Z') x p(z'|z)', (N by Z) x (Z by 1) = N by 1
        % Note: transpose ar_z_trans_condi from 1 by Z to Z by 1
        % Note: matrix multiply not dot multiply
        mt_evzp_condi_z = mt_val_cur * ar_z_trans_condi';

        % EVAL add on future utility, N by N + N by 1, broadcast again
        mt_utility = mt_utility + fl_beta*mt_evzp_condi_z;
        mt_utility(mt_c <= 0) = fl_nan_replace;

        % Optimization: remember matlab is column major, rows must be
        % choices, columns must be states
        % <https://en.wikipedia.org/wiki/Row-_and_column-major_order COLUMN-MAJOR>
        % mt_utility is N by N, rows are choices, cols are states.
        [ar_opti_val1_z, ar_opti_idx_z] = max(mt_utility);
        mt_val(:,it_z_i) = ar_opti_val1_z;
        mt_pol_a(:,it_z_i) = ar_a(ar_opti_idx_z);
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
% for choices outcomes, store as cell with two elements, first element is
% the y(a,z), outcome given states, the second element will be solved found
% in
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_ds_vf.html
% ff_ds_vf> and other distributions files. It stores what are the
% probability mass function of y, along with sorted unique values of y.

result_map = containers.Map('KeyType','char', 'ValueType','any');
result_map('mt_val') = mt_val;
result_map('mt_pol_idx') = mt_pol_idx;

result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
result_map('cl_mt_coh') = {f_coh(ar_z, ar_a'), zeros(1)};
result_map('cl_mt_pol_c') = {f_coh(ar_z, ar_a') - mt_pol_a, zeros(1)};
result_map('ar_st_pol_names') = ["mt_val", "cl_mt_pol_a", "cl_mt_coh", "cl_mt_pol_c"];

if (bl_post)
    bl_input_override = true;
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);
    result_map = ff_az_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
end

end
