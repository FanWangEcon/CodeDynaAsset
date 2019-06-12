%%
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

function result_map = ff_az_vf_vecsv(varargin)
%% FF_AZ_VF_VECSV solve infinite horizon exo shock + endo asset problem
% This program solves the infinite horizon dynamic single asset and single
% shock problem with vectorized codes.
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf.html
% ff_az_vf> shows looped codes.
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vec.html
% ff_az_vf_vec> shows vectorized codes. This file shows vectorized codes
% that is faster but is more memory intensive.
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
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/paramfunc/ffs_az_set_default_param.m ffs_az_set_default_param>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/paramfunc/ffs_az_get_funcgrid.m ffs_az_get_funcgrid>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/solvepost/ff_az_vf_post.m ff_az_vf_post>
%


%% Default
%
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button
%
% go to
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html
% ffs_az_set_default_param> to change parameters in param_map container.
% The parameters can also be updated here directly after obtaining them
% from ffs_az_set_default_param as we possibly change it_a_n and it_z_n
% here.

it_param_set = 4;
bl_input_override = true;
[param_map, support_map] = ffs_az_set_default_param(it_param_set);
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
st_func_name = 'ff_az_vf_vecsv';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, {'ar_a', 'mt_z_trans', 'ar_z'});
[ar_a, mt_z_trans, ar_z] = params_group{:};
% func_map
params_group = values(func_map, {'f_util_log', 'f_util_crra', 'f_cons'});
[f_util_log, f_util_crra, f_cons] = params_group{:};
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

% We did not need these in ff_az_vf or ff_az_vf_vec
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
    % Only this segment of code differs between ff_az_vf and ff_az_vf_vec
    % Store in cells results and retrieve, this is more memory intensive
    % than ff_az_vf_vec.

    % loop 1: over exogenous states
    for it_z_i = 1:length(ar_z)

        % Current Shock
        fl_z = ar_z(it_z_i);

        % Consumption and u(c) only need to be evaluated once
        if (it_iter == 1)

            % Consumption
            mt_c = f_cons(fl_z, ar_a, ar_a');

            % EVAL current utility: N by N, f_util defined earlier
            % slightly faster to explicitly write function
            if (fl_crra == 1)
                mt_utility = log(mt_c);
            else
                % slightly faster if write function here directly, but
                % speed gain is very small, more important to have single
                % location control of functions.
                mt_utility = f_util_crra(mt_c);
            end

            % Eliminate Complex Numbers
            mt_it_c_valid_idx = (mt_c <= 0);
            mt_utility(mt_it_c_valid_idx) = fl_nan_replace;

            % Store in cells
            cl_u_c_store{it_z_i} = mt_utility;
            cl_c_valid_idx{it_z_i} = mt_it_c_valid_idx;

        end

        % f(z'|z)
        ar_z_trans_condi = mt_z_trans(it_z_i,:);

        % EVAL EV((A',K'),Z'|Z) = V((A',K'),Z') x p(z'|z)', (N by Z) x (Z by 1) = N by 1
        mt_evzp_condi_z = mt_val_cur * ar_z_trans_condi';

        % EVAL add on future utility, N by N + N by 1
        mt_utility = cl_u_c_store{it_z_i} + fl_beta*mt_evzp_condi_z;

        % Index update
        % using the method below is much faster than index replace
        % see <https://fanwangecon.github.io/M4Econ/support/speed/index/fs_subscript.html fs_subscript>
        mt_it_c_valid_idx = cl_c_valid_idx{it_z_i};
        mt_utility = mt_utility.*(~mt_it_c_valid_idx) + fl_nan_replace*(mt_it_c_valid_idx);

        % Optimization: remember matlab is column major, rows must be
        % choices, columns must be states
        % <https://en.wikipedia.org/wiki/Row-_and_column-major_order COLUMN-MAJOR>
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

result_map = containers.Map('KeyType','char', 'ValueType','any');
result_map('mt_val') = mt_val;
result_map('mt_pol_a') = mt_pol_a;
result_map('mt_pol_idx') = mt_pol_idx;

if (bl_post)
    bl_input_override = true;
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);
    result_map = ff_az_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
end

end
