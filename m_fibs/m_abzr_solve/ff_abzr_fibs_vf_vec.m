%% Solve For+Inf+Borr+Save Dynamic Programming Problem (Vectorized)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function result_map = ff_abzr_fibs_vf_vec(varargin)
%% FF_ABZ_FIBS_VF_VEC solve infinite horizon exo shock + endo asset problem
% This program solves the infinite horizon dynamic single asset and single
% shock problem with vectorized codes.
% <https://fanwangecon.github.io/CodeDynaAsset/m_abzr/solve/html/ff_abzr_fibs_vf.html
% ff_abzr_fibs_vf> shows looped codes. The solution is the same.
%
% The model could be invoked mainly in sveral ways:
%
% # param_map('bl_default') = true;  param_map('bl_bridge') = false;
% param_map('bl_rollover') = true; Given these, default is possible, bridge
% loans are not needed because rollover is allowed for formal loans (or
% informal loans)
% # we change param_map('bl_bridge') = true, that means
% rollover is still allowed, but only allowed using informal sources,
% formal loans no longer allow for roll-over. Furthermore, if both
% bl_bridge and bl_rollover are false, that means we are not allowing for
% rollover at all, so households can not borrow such that they end up with
% negative cash-on-hand.
%
% Default simulation bl_bridge = false.
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
%    it_param_set = 4;
%    [param_map, support_map] = ffs_abzr_fibs_set_default_param(it_param_set);
%    % Chnage param_map keys for borrowing
%    param_map('fl_b_bd') = -20; % borrow bound
%    param_map('bl_default') = false; % true if allow for default
%    param_map('fl_c_min') = 0.0001; % u(c_min) when default
%    % Change Keys in param_map
%    param_map('it_a_n') = 75;
%    param_map('fl_z_r_borr_n') = 3;
%    param_map('it_z_wage_n') = 5;
%    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
%    param_map('fl_a_max') = 100;
%    param_map('fl_w') = 1.3;
%    % Change Keys support_map
%    support_map('bl_display') = false;
%    support_map('bl_post') = true;
%    support_map('bl_display_final') = false;
%    % Call Program with external parameters that override defaults.
%    ff_abzr_fibs_vf_vec(param_map, support_map);
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc/html/ffs_abzr_fibs_set_default_param.html ffs_abzr_fibs_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc/html/ffs_abzr_fibs_get_funcgrid.html ffs_abzr_fibs_get_funcgrid>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_min_c_cost_bridge.html ffs_fibs_min_c_cost_bridge>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_inf_bridge.html ffs_fibs_inf_bridge>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_min_c_cost.html ffs_fibs_min_c_cost>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_vf_post.html ff_az_vf_post>
%
% @seealso
%
% * for/inf + save + borr loop: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abzr_solve/html/ff_abzr_fibs_vf.html ff_abzr_fibs_vf>
% * for/inf + borr vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abzr_solve/html/ff_abzr_fibs_vf_vec.html ff_abzr_fibs_vf_vec>
% * for/inf + borr optimized-vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abzr_solve/html/ff_abzr_fibs_vf_vecsv.html ff_abzr_fibs_vf_vecsv>
%

%% Default
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button

it_param_set = 1;
bl_input_override = true;
[param_map, support_map] = ffs_abzr_fibs_set_default_param(it_param_set);

% Note: param_map and support_map can be adjusted here or outside to override defaults
% To generate results as if formal informal do not matter

% param_map('fl_r_fsv') = 0.025;
% param_map('fl_r_fbr') = 0.035;
% param_map('bl_b_is_principle') = false;
% param_map('st_forbrblk_type') = 'seg3';
% param_map('fl_forbrblk_brmost') = -19;
% param_map('fl_forbrblk_brleast') = -1;
% param_map('fl_forbrblk_gap') = -1.5;
% param_map('bl_b_is_principle') = false;
% param_map('it_a_n') = 750;
% param_map('fl_z_r_borr_n') = 5;
% param_map('it_z_wage_n') = 15;
% param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

[armt_map, func_map] = ffs_abzr_fibs_get_funcgrid(param_map, support_map); % 1 for override
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
    [armt_map, func_map] = ffs_abzr_fibs_get_funcgrid(param_map, support_map);
else
    % Override all
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

% append function name
st_func_name = 'ff_abzr_fibs_vf_vec';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, {'ar_a', 'mt_z_trans', 'ar_z_r_infbr_mesh_wage', 'ar_z_wage_mesh_r_infbr'});
[ar_a, mt_z_trans, ar_z_r_infbr_mesh_wage, ar_z_wage_mesh_r_infbr] = params_group{:};

% Formal choice Menu/Grid and Interest Rate Menu/Grid
params_group = values(armt_map, {'ar_forbrblk_r', 'ar_forbrblk'});
[ar_forbrblk_r, ar_forbrblk] = params_group{:};

% func_map
params_group = values(func_map, {'f_util_log', 'f_util_crra', 'f_coh', 'f_cons_coh_fbis', 'f_cons_coh_save'});
[f_util_log, f_util_crra, f_coh, f_cons_coh_fbis, f_cons_coh_save] = params_group{:};

% param_map
params_group = values(param_map, {'it_a_n', 'it_z_n', 'fl_crra', 'fl_beta', 'fl_c_min',...
    'fl_nan_replace', 'bl_default', 'bl_bridge', 'bl_rollover', 'fl_default_aprime'});
[it_a_n, it_z_n, fl_crra, fl_beta, fl_c_min, ...
    fl_nan_replace, bl_default, bl_bridge, bl_rollover, fl_default_aprime] = params_group{:};
params_group = values(param_map, {'it_maxiter_val', 'fl_tol_val', 'fl_tol_pol', 'it_tol_pol_nochange'});
[it_maxiter_val, fl_tol_val, fl_tol_pol, it_tol_pol_nochange] = params_group{:};

% param_map, Formal informal
params_group = values(param_map, {'fl_r_fsv', 'bl_b_is_principle'});
[fl_r_fsv, bl_b_is_principle] = params_group{:};

% support_map
params_group = values(support_map, {'bl_profile', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_display_minccost', 'bl_display_infbridge', ...
    'bl_time', 'bl_display_defparam', 'bl_display', 'it_display_every', 'bl_post'});
[bl_profile, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_display_minccost, bl_display_infbridge, ...
    bl_time, bl_display_defparam, bl_display, it_display_every, bl_post] = params_group{:};

%% Display Parameters

if (bl_display_defparam)
    fft_container_map_display(param_map);
    fft_container_map_display(support_map);
end

%% Initialize Output Matrixes
% include mt_pol_idx which we did not have in looped code

mt_val_cur = zeros(it_a_n,it_z_n);
mt_val = mt_val_cur - 1;
mt_pol_a = zeros(it_a_n,it_z_n);
mt_pol_a_cur = mt_pol_a - 1;
mt_pol_idx = zeros(it_a_n,it_z_n);
mt_pol_cons = zeros(it_a_n,it_z_n);

% collect optimal borrowing formal and informal choices
mt_pol_b_bridge = zeros(it_a_n,it_z_n);
mt_pol_inf_borr_nobridge = zeros(it_a_n,it_z_n);
mt_pol_for_borr = zeros(it_a_n,it_z_n);
mt_pol_for_save = zeros(it_a_n,it_z_n);

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

    % loop 1: over exogenous states
    for it_z_i = 1:it_z_n

        %% Solve the Formal Informal Problem for each a' and coh: c_forinf(a')
        % find the today's consumption maximizing formal and informal
        % choices given a' and coh. The formal and informal choices need to
        % generate exactly a', but depending on which formal and informal
        % joint choice is used, the consumption cost today a' is different.
        % Note here, a is principle + interests. Three areas:
        %
        % * *CASE A* a' > 0: savings, do not need to optimize over formal and
        % informal choices
        % * *CASE B* a' < 0 & coh < 0: need bridge loan to pay for unpaid debt, and
        % borrowing over-all, need to first pick bridge loan to pay for
        % debt, if bridge loan is insufficient, go into default. After
        % bridge loan, optimize over formal+informal, borrow+save joint
        % choices.
        % * *CASE C* a' < 0 & coh > 0: do not need to get informal bridge loans,
        % optimize over for+inf save, for+save+borr, inf+borr only, for
        % borrow only.
        %

        % 1. Current Shock
        fl_z_r_borr = ar_z_r_infbr_mesh_wage(it_z_i);
        fl_z_wage = ar_z_wage_mesh_r_infbr(it_z_i);

        % 2. cash-on-hand
        ar_coh = f_coh(fl_z_wage, ar_a);

        % 3. *CASE A* initiate consumption matrix as if all save
        mt_c = f_cons_coh_save(ar_coh, ar_a');

        % 3. if Bridge Loan is Needed

        % 4. *CASE B+C* get negative coh index and get borrowing choices index
        ar_coh_neg_idx = (ar_coh <= 0);

        ar_a_neg_idx = (ar_a < 0);
        ar_coh_neg = ar_coh(ar_coh_neg_idx);
        ar_a_neg = ar_a(ar_a_neg_idx);

        % 5. if coh > 0 and ap < 0, can allow same for+inf result to all coh.
        % The procedure below works regardless of how ar_coh is sorted. get
        % the index of all negative coh elements as well as first
        % non-negative element. We solve the formal and informal problem at
        % these points, note that we only need to solve the formal and
        % informal problem for positive coh level once.
        ar_coh_first_pos_idx = (cumsum(ar_coh_neg_idx == 0) == 1);
        ar_coh_forinfsolve_idx = (ar_coh_first_pos_idx | ar_coh_neg_idx);
        ar_coh_forinfsolve_a_neg_idx = (ar_coh(ar_coh_forinfsolve_idx) <= 0);

        % 6. *CASE B + C* Negative asset choices (borrowing), 1 col Case C
        % negp1: negative coh + 1, 1 meaning 1 positive coh, first positive
        % coh column index element grabbed.
        mt_coh_negp1_mesh_neg_aprime = zeros(size(ar_a_neg')) + ar_coh(ar_coh_forinfsolve_idx);
        mt_neg_aprime_mesh_coh_negp1 = zeros(size(mt_coh_negp1_mesh_neg_aprime)) + ar_a_neg';

        if (bl_bridge)
            %         mt_neg_aprime_mesh_coh_1col4poscoh = zeros([length(ar_a_neg), (length(ar_coh_neg)+1)]) + ar_a_neg';
            %         ar_coh_neg_idx_1col4poscoh = ar_coh_neg_idx(1:(length(ar_coh_neg)+1));

            % 6. *CASE B* Solve for: if (fl_ap < 0) and if (fl_coh < 0)
            [mt_aprime_nobridge_negcoh, ~, mt_c_bridge_negcoh] = ffs_fibs_inf_bridge(...
                bl_b_is_principle, fl_z_r_borr, ...
                mt_neg_aprime_mesh_coh_negp1(:,ar_coh_forinfsolve_a_neg_idx), ...
                mt_coh_negp1_mesh_neg_aprime(:,ar_coh_forinfsolve_a_neg_idx), ...
                bl_display_infbridge, bl_input_override);

            % generate mt_aprime_nobridge
            mt_neg_aprime_mesh_coh_negp1(:, ar_coh_forinfsolve_a_neg_idx) = mt_aprime_nobridge_negcoh;
        else
            % no bridge loan needed means roll over is allowed.
            mt_neg_aprime_mesh_coh_negp1 = ar_a_neg';
        end

        % 7. *CASE B + C* formal and informal joint choices, 1 col Case C
        bl_input_override = true;
        [ar_max_c_nobridge, ~, ~, ~] = ...
            ffs_fibs_min_c_cost(...
            bl_b_is_principle, fl_z_r_borr, fl_r_fsv, ...
            ar_forbrblk_r, ar_forbrblk, ...
            mt_neg_aprime_mesh_coh_negp1(:), ...
            bl_display_minccost, bl_input_override);

        %% Update Consumption Matrix *CASE A + B + C* Consumptions
        % Current mt_c is assuming all to be case A
        %
        % * Update Columns for case B (negative coh)
        % * Update Columns for case C (1 column): ar_coh_first_pos_idx,
        % included in ar_coh_forinfsolve_idx
        % * Update Columns for all case C: ~ar_coh_neg_idx using 1 column
        % result
        %

        % 1. Initalize all Neg Aprime consumption cost of aprime inputs
        % Initialize
        mt_max_c_nobridge_a_neg = zeros([length(ar_a_neg), length(ar_coh)]) + 0;
        mt_c_bridge_coh_a_neg = zeros(size(mt_max_c_nobridge_a_neg)) + 0;

        % 2. Fill in *Case B* and *Case C* (one column) Other C-cost
        mt_max_c_nobridge_negcohp1 = reshape(ar_max_c_nobridge, [size(mt_neg_aprime_mesh_coh_negp1)]);

        if (bl_bridge)
            % 2. Fill in *Case B* Bridge C-cost
            mt_c_bridge_coh_a_neg(:, ar_coh_neg_idx) = mt_c_bridge_negcoh;

            % 2. Fill in *Case B* and *Case C* (one column) Other C-cost
            mt_max_c_nobridge_a_neg(:, ar_coh_forinfsolve_idx) = mt_max_c_nobridge_negcohp1;
            mt_max_c_nobridge_a_neg(:, ~ar_coh_forinfsolve_idx) = ...
                zeros(size(mt_c(ar_a_neg_idx, ~ar_coh_forinfsolve_idx))) ...
                + mt_max_c_nobridge_negcohp1(:, ~ar_coh_forinfsolve_a_neg_idx);
        else
            mt_max_c_nobridge_a_neg = zeros([length(ar_a_neg), length(ar_coh)]) + mt_max_c_nobridge_negcohp1;
        end

        % 3. Consumption for B + C Cases
        % note, the c cost of aprime is the same for all coh > 0, but mt_c
        % is different still for each coh and aprime.
        mt_c_forinfsolve = f_cons_coh_fbis(ar_coh, mt_c_bridge_coh_a_neg + mt_max_c_nobridge_a_neg);

        % 4. Update with Case B and C
        mt_c(ar_a_neg_idx, :) = mt_c_forinfsolve;

        %% Solve Optimization Problem: max_{a'} (u(c_forinf(a')) + EV(a',z'))
        % 1. EVAL current utility: N by N, f_util defined earlier
        if (fl_crra == 1)
            mt_utility = f_util_log(mt_c);
            fl_u_cmin = f_util_log(fl_c_min);
        else
            mt_utility = f_util_crra(mt_c);
            fl_u_cmin = f_util_crra(fl_c_min);
        end

        % 2. f(z'|z)
        ar_z_trans_condi = mt_z_trans(it_z_i,:);

        % 3. EVAL EV((A',K'),Z'|Z) = V((A',K'),Z') x p(z'|z)', (N by Z) x (Z by 1) = N by 1
        % Note: transpose ar_z_trans_condi from 1 by Z to Z by 1
        % Note: matrix multiply not dot multiply
        mt_evzp_condi_z = mt_val_cur * ar_z_trans_condi';

        % 4. EVAL add on future utility, N by N + N by 1, broadcast again
        mt_utility = mt_utility + fl_beta*mt_evzp_condi_z;

        if (bl_default)
            % if default: only today u(cmin), transition out next period, debt wiped out
            mt_utility(mt_c <= fl_c_min) = fl_u_cmin + fl_beta*mt_evzp_condi_z(ar_a == fl_default_aprime);
        else
            % if default is not allowed: v = fl_nan_replace
            mt_utility(mt_c <= fl_c_min) = fl_nan_replace;
        end

        % Set below threshold c to c_min
        mt_c(mt_c < fl_c_min) = fl_c_min;

        % 5. no bridge and no rollover allowed
        if( ~bl_rollover && ~bl_bridge)
            if (bl_default)
                % if default: only today u(cmin), transition out next period, debt wiped out
                mt_utility(:, ar_coh_neg_idx) = fl_u_cmin + fl_beta*mt_evzp_condi_z(ar_a == fl_default_aprime);
            else
                % if default is not allowed: v = fl_nan_replace
                mt_utility(:, ar_coh_neg_idx) = fl_nan_replace;
            end
        end

        % 5. Optimization: remember matlab is column major, rows must be
        % choices, columns must be states
        % <https://en.wikipedia.org/wiki/Row-_and_column-major_order COLUMN-MAJOR>
        % mt_utility is N by N, rows are choices, cols are states.
        [ar_opti_val_z, ar_opti_idx_z] = max(mt_utility);
        [it_choies_n, it_states_n] = size(mt_utility);
        ar_add_grid = linspace(0, it_choies_n*(it_states_n-1), it_states_n);
        ar_opti_linear_idx_z = ar_opti_idx_z + ar_add_grid;
        ar_opti_aprime_z = ar_a(ar_opti_idx_z);
        ar_opti_c_z = mt_c(ar_opti_linear_idx_z);

        % 6. Handle Default is optimal or not
        if (bl_default)
            % if defaulting is optimal choice, at these states, not required
            % to default, non-default possible, but default could be optimal
            ar_opti_aprime_z(ar_opti_c_z <= fl_c_min) = fl_default_aprime;
            ar_opti_idx_z(ar_opti_c_z <= fl_c_min) = find(ar_a == fl_default_aprime);
        else
            % if default is not allowed, then next period same state as now
            % this is absorbing state, this is the limiting case, single
            % state space point, lowest a and lowest shock has this.
            ar_opti_aprime_z(ar_opti_c_z <= fl_c_min) = min(ar_a);
        end

        % 6. no bridge and no rollover allowed
        if( ~bl_rollover && ~bl_bridge)
            if (bl_default)
                % if default: only today u(cmin), transition out next period, debt wiped out
                ar_opti_aprime_z(ar_coh_neg_idx) = fl_default_aprime;
            else
                % if default is not allowed: v = fl_nan_replace
                ar_opti_aprime_z(ar_coh_neg_idx) = ar_a(fl_nan_replace);
            end
        end

        %% Store Optimal Choices Current Iteration
        mt_val(:,it_z_i) = ar_opti_val_z;
        mt_pol_a(:,it_z_i) = ar_opti_aprime_z;
        mt_pol_cons(:,it_z_i) = ar_opti_c_z;

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

% Find optimal Formal Informal Choices. Could have saved earlier, but was
% wasteful of resources
for it_z_i = 1:it_z_n
    for it_a_j = 1:it_a_n
        fl_z_r_borr = ar_z_r_infbr_mesh_wage(it_z_i);
        fl_z_wage = ar_z_wage_mesh_r_infbr(it_z_i);        

        fl_a = ar_a(it_a_j);
        fl_coh = f_coh(fl_z_wage, fl_a);
        fl_a_opti = mt_pol_a(it_a_j, it_z_i);

        param_map('fl_r_inf') = fl_z_r_borr;

        % call formal and informal function.
        [~, fl_opti_b_bridge, fl_opti_inf_borr_nobridge, fl_opti_for_borr, fl_opti_for_save] = ...
            ffs_fibs_min_c_cost_bridge(fl_a_opti, fl_coh, ...
            param_map, support_map, armt_map, func_map, bl_input_override);

        % store savings and borrowing formal and inf optimal choices
        mt_pol_b_bridge(it_a_j,it_z_i) = fl_opti_b_bridge;
        mt_pol_inf_borr_nobridge(it_a_j,it_z_i) = fl_opti_inf_borr_nobridge;
        mt_pol_for_borr(it_a_j,it_z_i) = fl_opti_for_borr;
        mt_pol_for_save(it_a_j,it_z_i) = fl_opti_for_save;

    end
end

result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
result_map('cl_mt_coh') = {f_coh(ar_z_r_infbr_mesh_wage, ar_a'), zeros(1)};

result_map('cl_mt_pol_c') = {mt_pol_cons, zeros(1)};
result_map('cl_mt_pol_b_bridge') = {mt_pol_b_bridge, zeros(1)};
result_map('cl_mt_pol_inf_borr_nobridge') = {mt_pol_inf_borr_nobridge, zeros(1)};
result_map('cl_mt_pol_for_borr') = {mt_pol_for_borr, zeros(1)};
result_map('cl_mt_pol_for_save') = {mt_pol_for_save, zeros(1)};

result_map('ar_st_pol_names') = ["cl_mt_pol_a", "cl_mt_pol_coh", "cl_mt_pol_c", ...
    "cl_mt_pol_b_bridge", "cl_mt_pol_inf_borr_nobridge", "cl_mt_pol_for_borr", "cl_mt_pol_for_save"];

% Get Discrete Choice Outcomes
result_map = ffs_fibs_identify_discrete(result_map, bl_input_override);

%% Post Solution Graph and Table Generation
% Note in comparison with *abzr*, results here, even when using identical
% parameters would differ because in *abzr* solved where choices are
% principle. Here choices are principle + interests in order to facilitate
% using the informal choice functions.
%
% Note that this means two things are
% different, on the one hand, the value of asset for to coh is different
% based on the grid of assets. If the asset grid is negative, now per grid
% point, there is more coh because that grid point of asset no longer has
% interest rates. On the other hand, if one has positive asset grid point
% on arrival, that is worth less to coh. Additionally, when making choices
% for the next period, now choices aprime includes interests. What these
% mean is that the a grid no longer has the same meaning. We should expect
% at higher savings levels, for the same grid points, if optimal grid
% choices are the same as before, consumption should be lower when b
% includes interest rates and principle. This is however, not true when
% arriving in a period with negative a levels, for the same negative a
% level and same a prime negative choice, could have higher consumption
% here becasue have to pay less interests on debt. This tends to happen for
% smaller levels of borrowing choices.
%
% Graphically, when using interest + principle, big difference in
% consumption as a fraction of (coh - aprime) figure. In those figures,
% when counting in principles only, the gap in coh and aprime is
% consumption, but now, as more is borrowed only a small fraction of coh
% and aprime gap is consumption, becuase aprime/(1+r) is put into
% consumption.

if (bl_post)
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);

    % Standard AZ graphs
    result_map = ff_az_vf_post(param_map, support_map, armt_map, func_map, result_map);

    % Graphs for results_map with FIBS contents
    result_map = ff_az_fibs_vf_post(param_map, support_map, armt_map, func_map, result_map);
end

end
