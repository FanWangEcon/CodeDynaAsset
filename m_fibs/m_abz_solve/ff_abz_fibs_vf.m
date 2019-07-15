%% Solve For+Inf+Borr+Save Dynamic Programming Problem (Loop)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function result_map = ff_abz_fibs_vf(varargin)
%% FF_ABZR_FIBS_VF borr + save one asset formal informal + loop
% This program solves the infinite horizon dynamic single asset and single
% shock problem with loops. This file contains codes that processes
% borrowing and handles formal and informal choices. R shock.
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
% * mt_cons matrix states_n by shock_n matrix of optimal consumption
% levels, unlike modele without formal and informal choices, where we know
% c from coh and a, here this needed to be stored because it is the results
% from with joint category maximization problem.
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
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_paramfunc/html/ffs_abz_fibs_set_default_param.html ffs_abz_fibs_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_paramfunc/html/ffs_abz_fibs_get_funcgrid.html ffs_abz_fibs_get_funcgrid>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_min_c_cost_bridge.html ffs_fibs_min_c_cost_bridge>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_inf_bridge.html ffs_fibs_inf_bridge>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_min_c_cost.html ffs_fibs_min_c_cost>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_vf_post.html ff_az_vf_post>
%
% @seealso
%
% * for/inf + save + borr loop: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_solve/html/ff_abz_fibs_vf.html ff_abz_fibs_vf>
% * for/inf + borr vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_solve/html/ff_abz_fibs_vf_vec.html ff_abz_fibs_vf_vec>
% * for/inf + borr optimized-vectorized: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_solve/html/ff_abz_fibs_vf_vecsv.html ff_abz_fibs_vf_vecsv>
%

%% Default
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button

it_param_set = 4;
bl_input_override = true;
[param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);

% Note: param_map and support_map can be adjusted here or outside to override defaults
% To generate results as if formal informal do not matter
param_map('it_a_n') = 35;
param_map('it_z_n') = 7;
param_map('it_maxiter_val') = 20;
% param_map('fl_r_fsv') = 0.025;
% param_map('fl_r_inf') = 0.035;
% param_map('fl_r_inf_bridge') = 0.035;
% param_map('fl_r_fbr') = 0.035;
% param_map('bl_b_is_principle') = false;
% param_map('st_forbrblk_type') = 'seg3';
% param_map('fl_forbrblk_brmost') = -19;
% param_map('fl_forbrblk_brleast') = -1;
% param_map('fl_forbrblk_gap') = -1.5;
% param_map('bl_b_is_principle') = false;

[armt_map, func_map] = ffs_abz_fibs_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
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
    [armt_map, func_map] = ffs_abz_fibs_get_funcgrid(param_map, support_map, bl_input_override);
else
    % Override all
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

% append function name
st_func_name = 'ff_abz_fibs_vf';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, {'ar_a', 'mt_z_trans', 'ar_z'});
[ar_a, mt_z_trans, ar_z] = params_group{:};

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
params_group = values(param_map, {'fl_r_inf', 'fl_r_fsv', 'bl_b_is_principle'});
[fl_r_inf, fl_r_fsv, bl_b_is_principle] = params_group{:};

% support_map
params_group = values(support_map, {'bl_profile', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_display_minccost', 'bl_display_infbridge', ...
    'bl_time', 'bl_display', 'it_display_every', 'bl_post'});
[bl_profile, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_display_minccost, bl_display_infbridge, ...
    bl_time, bl_display, it_display_every, bl_post] = params_group{:};

%% Initialize Output Matrixes

mt_val_cur = zeros(length(ar_a),length(ar_z));
mt_val = mt_val_cur - 1;
mt_pol_a = zeros(length(ar_a),length(ar_z));
mt_pol_a_cur = mt_pol_a - 1;
mt_pol_cons = zeros(length(ar_a),length(ar_z));

% collect optimal borrowing formal and informal choices
mt_pol_b_bridge = zeros(length(ar_a),length(ar_z));
mt_pol_inf_borr_nobridge = zeros(length(ar_a),length(ar_z));
mt_pol_for_borr = zeros(length(ar_a),length(ar_z));
mt_pol_for_save = zeros(length(ar_a),length(ar_z));

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
    % handling borrowing and default possibility

    % loop 1: over exogenous states
    for it_z_i = 1:length(ar_z)
        fl_z = ar_z(it_z_i);

        % loop 2: over endogenous states
        for it_a_j = 1:length(ar_a)

            % Get asset state
            fl_a = ar_a(it_a_j);

            % Initialize storage
            ar_val_cur = zeros(size(ar_a));
            ar_c_cur = zeros(size(ar_a));
            ar_b_bridge = zeros(size(ar_a));
            ar_inf_borr_nobridge = zeros(size(ar_a));
            ar_for_borr = zeros(size(ar_a));
            ar_for_save = zeros(size(ar_a));

            % calculate cash on hand
            fl_coh = f_coh(fl_z, fl_a);

            % loop 3: over choices
            for it_ap_k = 1:length(ar_a)

                % get next period asset choice
                fl_ap = ar_a(it_ap_k);

                %% Compute Consumption given Borrowing and Savings
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
                % * *CASE C* $ a' < 0 & coh > 0: do not need to get informal bridge loans,
                % optimize over for+inf save, for+save+borr, inf+borr only, for
                % borrow only.
                %

                if (fl_ap < 0)

                    % Calculate Bridge Loan Borrowing
                    if (bl_bridge && fl_coh < 0)

                        bl_input_override = true;
                        [fl_aprime_nobridge, fl_b_bridge, fl_c_bridge] = ffs_fibs_inf_bridge(...
                            bl_b_is_principle, fl_r_inf, fl_ap, fl_coh, ...
                            bl_display_infbridge, bl_input_override);

                    else

                        fl_aprime_nobridge = fl_ap;
                        fl_b_bridge = 0;
                        fl_c_bridge = 0;

                    end

                    % Find Optimal Formal Informal Borrow Save Combo
                    % calculate consumption gain from formal + informal
                    % borrowing and savings choices.
                    bl_input_override = true;
                    [fl_max_c_nobridge, fl_inf_borr_nobridge, fl_for_borr, fl_for_save] = ...
                        ffs_fibs_min_c_cost(...
                        bl_b_is_principle, fl_r_inf, fl_r_fsv, ...
                        ar_forbrblk_r, ar_forbrblk, ...
                        fl_aprime_nobridge, bl_display_minccost, bl_input_override);

                    % Compute Consumption given Formal and Informal joint
                    % consumption with formal borrow menu + bridge loans.
                    fl_c = f_cons_coh_fbis(fl_coh, fl_max_c_nobridge + fl_c_bridge);

                else

                    % consumption with savings
                    fl_c = f_cons_coh_save(fl_coh, fl_ap);

                    % assign values for formal and informal choices
                    % possible that fl_coh < 0, but if then fl_ap > 0 is
                    % not valid choice
                    [fl_b_bridge, fl_inf_borr_nobridge, fl_for_borr, fl_for_save] = deal(0, 0, 0, fl_ap);

                end

                %% Compute Utility With Default
                % if rollover is not allowed and bridge is not allowed,
                % then as long as coh <= 0, also treat as not allowed
                % states.
                % assign u(c)
                if (fl_c <= fl_c_min || ...
                    ( ~bl_rollover && ~bl_bridge && fl_coh < fl_c_min))

                    if (bl_default)
                        % defaults
                        % current utility: only today u(cmin)
                        ar_val_cur(it_ap_k) = fl_u_cmin;
                        % transition out next period, debt wiped out
                        for it_az_q = 1:length(ar_z)
                            ar_val_cur(it_ap_k) = ar_val_cur(it_ap_k) + ...
                                fl_beta*mt_z_trans(it_z_i, it_az_q)*mt_val_cur((ar_a == fl_default_aprime), it_az_q);
                        end

                        % Replace Consumption if default cmin
                        fl_c = fl_c_min;
                    else
                        % if default is not allowed: v = fl_nan_replace
                        ar_val_cur(it_ap_k) = fl_nan_replace;

                        % Replace Consumption if no default nan
                        fl_c = 0;
                    end

                    % no action, defaulting
                    fl_b_bridge = 0;
                    fl_inf_borr_nobridge = 0;
                    fl_for_borr = 0;
                    fl_for_save = 0;

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

                %% Store Values

                % Could get the formal and informal values from
                % ffs_fibs_min_c_cost_bridge.m
%                 bl_input_override = true;
%                 [fl_c, fl_b_bridge, fl_inf_borr_nobridge, fl_for_borr, fl_for_save] = ...
%                     ffs_fibs_min_c_cost_bridge(fl_ap, fl_coh, ...
%                     param_map, support_map, armt_map, func_map, bl_input_override);

                % Store consumption
                ar_c_cur(it_ap_k) = fl_c;

                % Save/Update Borrowing Information
                ar_b_bridge(it_ap_k) = fl_b_bridge;
                ar_inf_borr_nobridge(it_ap_k) = fl_inf_borr_nobridge;
                ar_for_borr(it_ap_k) = fl_for_borr;
                ar_for_save(it_ap_k) = fl_for_save;

            end

            %% Optimize over Next Period Asset Choices
            % optimal choice value
            [fl_opti_val_z, fl_opti_idx_z] = max(ar_val_cur);
            fl_opti_aprime_z = ar_a(fl_opti_idx_z);
            fl_opti_c_z = ar_c_cur(fl_opti_idx_z);

            % corresponding optimal borrowing and savings choices
            fl_opti_b_bridge = ar_b_bridge(fl_opti_idx_z);
            fl_opti_inf_borr_nobridge = ar_inf_borr_nobridge(fl_opti_idx_z);
            fl_opti_for_borr = ar_for_borr(fl_opti_idx_z);
            fl_opti_for_save = ar_for_save(fl_opti_idx_z);

            %% Find Optimal Choices for Defaults or Not
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
                    fl_opti_aprime_z = min(ar_a);
                end
            end

            %% Store Optimal Choices and Value Given(a,z)

            % store overal savings, value and consumption
            mt_val(it_a_j,it_z_i) = fl_opti_val_z;
            mt_pol_a(it_a_j,it_z_i) = fl_opti_aprime_z;
            mt_pol_cons(it_a_j,it_z_i) = fl_opti_c_z;

            % store savings and borrowing formal and inf optimal choices
            mt_pol_b_bridge(it_a_j,it_z_i) = fl_opti_b_bridge;
            mt_pol_inf_borr_nobridge(it_a_j,it_z_i) = fl_opti_inf_borr_nobridge;
            mt_pol_for_borr(it_a_j,it_z_i) = fl_opti_for_borr;
            mt_pol_for_save(it_a_j,it_z_i) = fl_opti_for_save;

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

result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
result_map('cl_mt_coh') = {f_coh(ar_z, ar_a'), zeros(1)};

result_map('cl_mt_pol_c') = {mt_pol_cons, zeros(1)};
result_map('cl_mt_pol_b_bridge') = {mt_pol_b_bridge, zeros(1)};
result_map('cl_mt_pol_inf_borr_nobridge') = {mt_pol_inf_borr_nobridge, zeros(1)};
result_map('cl_mt_pol_for_borr') = {mt_pol_for_borr, zeros(1)};
result_map('cl_mt_pol_for_save') = {mt_pol_for_save, zeros(1)};

result_map('ar_st_pol_names') = ["cl_mt_pol_a", "cl_mt_coh", "cl_mt_pol_c", ...
    "cl_mt_pol_b_bridge", "cl_mt_pol_inf_borr_nobridge", "cl_mt_pol_for_borr", "cl_mt_pol_for_save"];

% Get Discrete Choice Outcomes
result_map = ffs_fibs_identify_discrete(result_map, bl_input_override);

%% Post Solution Graph and Table Generation
% Note in comparison with *abz*, results here, even when using identical
% parameters would differ because in *abz* solved where choices are
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
    bl_input_override = true;
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);

    % Standard AZ graphs
    result_map = ff_az_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Graphs for results_map with FIBS contents
    result_map = ff_az_fibs_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
end

%% Display Various Containers
bl_display_defparam = true;
if (bl_display_defparam)
    
    %% Display 1 support_map    
    fft_container_map_display(support_map);
        
    %% Display 2 armt_map
    fft_container_map_display(armt_map);

    %% Display 3 param_map
    fft_container_map_display(param_map);
    
    %% Display 4 func_map
    fft_container_map_display(func_map);
    
    %% Display 5 result_map
    fft_container_map_display(result_map);
    
end

end
