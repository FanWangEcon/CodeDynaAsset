%% Solve Save + Borr Dynamic Programming Problem (Loop)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function result_map = ff_abz_vf(varargin)
%% FF_ABZ_VF solve infinite horizon exo shock + endo asset problem
% This program solves the infinite horizon dynamic single asset and single
% shock problem with loops. This file contains codes that processes
% borrowing.
%
% The borrowing problem is very similar to the savings problem. The code
% could be identical if one does not have to deal with default. The main
% addition here in comparison to the savings only code
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf.html
% ff_az_vf> is the ability to deal with default. 
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
%    param_map('it_z_n') = 11;
%    param_map('fl_a_max') = 100;
%    param_map('fl_w') = 1.3;
%    % Change Keys support_map
%    support_map('bl_display') = false;
%    support_map('bl_post') = true;
%    support_map('bl_display_final') = false;
%    % Call Program with external parameters that override defaults.
%    ff_abz_vf(param_map, support_map);
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
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button

it_param_set = 1;
bl_input_override = true;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% Note: param_map and support_map can be adjusted here or outside to override defaults
% param_map('it_a_n') = 750;
% param_map('it_z_n') = 15;
% param_map('fl_r_save') = 0.025;
% param_map('fl_r_borr') = 0.035;

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

% append function name
st_func_name = 'ff_abz_vf';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, {'ar_a', 'mt_z_trans', 'ar_z'});
[ar_a, mt_z_trans, ar_z] = params_group{:};

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
    'bl_time', 'bl_display', 'it_display_every', 'bl_post'});
[bl_profile, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display, it_display_every, bl_post] = params_group{:};

%% Initialize Output Matrixes

mt_val_cur = zeros(length(ar_a),length(ar_z));
mt_val = mt_val_cur - 1;
mt_pol_a = zeros(length(ar_a),length(ar_z));
mt_pol_a_cur = mt_pol_a - 1;

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

% Utility at-Default/at-limiting-case-when-nodefault
if (fl_crra == 1)
    fl_u_cmin = f_util_log(fl_c_min);
else
    fl_u_cmin = f_util_crra(fl_c_min);
end


% Value Function Iteration
while bl_vfi_continue
    it_iter = it_iter + 1;

    %% Iterate over a and z states
    % Solve Problem if Default: coh < borr_bound, then: u(c_min) + EV(a'=0,z')
    % these are the values for defaulter
    % if do not borrow enough to pay debt: a'=0
    % if save more than what you have: a'=0
    % if borrowing, possible that all ar_val_cur values are
    % u(cmin), when that is the case, utility equal to
    % c_min, this means given current choice set, can only default,
    % get utility at c_min level, and a' go to 0. See
    % maximization over loop 3 choices for loop 1+2 states. see
    % <https://fanwangecon.github.io/CodeDynaAsset/docs/README_cminymin_borrsave.html
    % README_cminymin_borrsave> for additional discussions.
    %
    % when default is not allowed, there should be exactly one
    % element of the state space where we enter this condition.
    % That is at a_min and z_min. At this point, the natural
    % borrowing constraint when fully used up leads to c = 0.
    % When default is not allowed, this is the limiting case.
    % We solve the limiting case as if default is possible. see
    % the figure here:
    % <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ffs_abz_get_funcgrid/test_borr/html/ffs_abz_get_funcgrid_nodefault.html
    % ffs_abz_get_funcgrid_nodefault>. For the figure under
    % section *Generate Borrowing A Grid with Default*, there are
    % no points to the lower right of the red horizontal and
    % verticle lines, because no default is allowed. But the
    % intersection of the two red lines is this limiting case
    % where c = 0.
    %
    % when default is allowed. See again
    % <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ffs_abz_get_funcgrid/test_borr/html/ffs_abz_get_funcgrid_nodefault.html
    % ffs_abz_get_funcgrid_nodefault>. Again under section
    % *Generate Borrowing A Grid with Default*, see that there
    % is an area to the right bottom of the two blue horizontal
    % and vertical lines. For all points in that area, the
    % household has to default, they can not borrow even at max
    % to get positive consumption. At other points higher than
    % the blue horizontal line, it is also possible for
    % defaulting to be the optimal choice.

    % loop 1: over exogenous states
    for it_z_i = 1:length(ar_z)
        fl_z = ar_z(it_z_i);

        % loop 2: over endogenous states
        for it_a_j = 1:length(ar_a)
            fl_a = ar_a(it_a_j);
            ar_val_cur = zeros(size(ar_a));

            % calculate cash on hand
            fl_coh = f_coh(fl_z, fl_a);

            % loop 3: over choices
            for it_ap_k = 1:length(ar_a)

                % get next period asset choice
                fl_ap = ar_a(it_ap_k);

                % calculate consumption
                fl_c = f_cons_coh(fl_coh, fl_ap);

                % assign u(c)
                if (fl_c <= fl_c_min)
                    if (bl_default)
                        % defaults
                        % current utility: only today u(cmin)
                        ar_val_cur(it_ap_k) = fl_u_cmin;
                        % transition out next period, debt wiped out
                        for it_az_q = 1:length(ar_z)
                            ar_val_cur(it_ap_k) = ar_val_cur(it_ap_k) + ...
                                fl_beta*mt_z_trans(it_z_i, it_az_q)*mt_val_cur((ar_a == fl_default_aprime), it_az_q);
                        end
                    else
                        % if default is not allowed: v = fl_nan_replace
                        ar_val_cur(it_ap_k) = fl_nan_replace;
                    end
                else
                    % Solve Optimization Problem: max_{a'} u(c(a,a',z)) + beta*EV(a',z')
                    % borrowed enough to pay debt (and borrowing limit not exceeded)
                    % saved only the coh available.
                    % current utility
                    if (fl_crra == 1)
                        ar_val_cur(it_ap_k) = f_util_log(fl_c);
                    else
                        ar_val_cur(it_ap_k) = f_util_crra(fl_c);
                    end
                    % loop 4: add future utility, integration--loop over future shocks
                    for it_az_q = 1:length(ar_z)
                        ar_val_cur(it_ap_k) = ar_val_cur(it_ap_k) + ...
                            fl_beta*mt_z_trans(it_z_i, it_az_q)*mt_val_cur(it_ap_k, it_az_q);
                    end
                end
            end

            % optimal choice value
            [fl_opti_val_z, fl_opti_idx_z] = max(ar_val_cur);
            fl_opti_aprime_z = ar_a(fl_opti_idx_z);
            fl_opti_c_z = f_cons_coh(fl_coh, fl_opti_aprime_z);

            % Handle Default is optimal or not
            if (fl_opti_c_z <= fl_c_min)
                if (bl_default)
                    % if defaulting is optimal choice, at these states, not required
                    % to default, non-default possible, but default could be optimal
                    fl_opti_aprime_z = fl_default_aprime;
                else
                    % if default is not allowed, then next period same state as now
                    % this is absorbing state, this is the limiting case, single
                    % state space point, lowest a and lowest shock has this.
                    fl_opti_aprime_z = fl_a;
                end
            end

            %% Store Optimal Choices and Value Given(a,z)
            mt_val(it_a_j,it_z_i) = fl_opti_val_z;
            mt_pol_a(it_a_j,it_z_i) = fl_opti_aprime_z;

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

result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
result_map('cl_mt_pol_coh') = {f_coh(ar_z, ar_a'), zeros(1)};
result_map('cl_mt_pol_c') = {f_cons_checkcmin(ar_z, ar_a', mt_pol_a), zeros(1)};
result_map('ar_st_pol_names') = ["cl_mt_pol_a", "cl_mt_pol_coh", "cl_mt_pol_c"];

if (bl_post)
    bl_input_override = true;
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);
    result_map = ff_az_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
end

end
