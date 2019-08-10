%% Risky + Safe Asset (Save + Borr + FIBS) Interpolated-Percentage (Optimized-Vectorized)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function result_map = ff_ipwkbz_fibs_vf_vecsv(varargin)
%% FF_IPWKBZ_VF_VECSV solve infinite horizon exo shock + endo asset problem
% This is a modified version of
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbz/solve/html/ff_ipwkbz_vf_vecsv.html
% ff_ipwkbz_vf_vecsv>, to see how this function solves the formal and
% savings risky and safe asset problem with formal and informal choices,
% compare the code here and from
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbz/solve/html/ff_ipwkbz_vf_vecsv.html
% ff_ipwkbz_vf_vecsv> side by side.
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
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkbz/paramfunc/ff_ipwkbz_evf.m ff_ipwkbz_evf>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkbz/paramfunc/ffs_ipwkbz_set_default_param.m ffs_ipwkbz_set_default_param>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_ipwkbz/paramfunc/ffs_ipwkbz_get_funcgrid.m ffs_ipwkbz_get_funcgrid>
% * <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solvepost/ff_akz_vf_post.m ff_akz_vf_post>
%

%% Default
% * it_param_set = 1: quick test
% * it_param_set = 2: benchmark run
% * it_param_set = 3: benchmark profile
% * it_param_set = 4: press publish button

it_param_set = 4;
[param_map, support_map] = ffs_ipwkbz_fibs_set_default_param(it_param_set);

% parameters can be set inside ffs_ipwkbz_set_default_param or updated here
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
params_len = length(varargin);
if params_len <= 2
    [armt_map, func_map] = ffs_ipwkbz_fibs_get_funcgrid(param_map, support_map); % 1 for override
    default_params = {param_map support_map armt_map func_map};
end

%% Parse Parameters 1

% if varargin only has param_map and support_map,
[default_params{1:params_len}] = varargin{:};
param_map = [param_map; default_params{1}];
support_map = [support_map; default_params{2}];
if params_len >= 1 && params_len <= 2
    % If override param_map, re-generate armt and func if they are not
    % provided
    [armt_map, func_map] = ffs_ipwkbz_fibs_get_funcgrid(param_map, support_map);
else
    % Override all
    armt_map = default_params{3};
    func_map = default_params{4};
end

% append function name
st_func_name = 'ff_ipwkbz_fibs_vf_vecsv';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Parse Parameters 2

% armt_map
params_group = values(armt_map, ...
    {'ar_w_perc', 'ar_w_level_full', 'ar_coh_bridge_perc', 'ar_z'});
[ar_w_perc, ar_w_level_full, ar_coh_bridge_perc, ar_z] = params_group{:};
params_group = values(armt_map, {'ar_interp_c_grid', 'ar_interp_coh_grid', ...
    'ar_a_meshk', 'ar_k_mesha', ...
    'mt_interp_coh_grid_mesh_z', 'mt_z_mesh_coh_interp_grid',...
    'mt_interp_coh_grid_mesh_w_perc',...
    'mt_w_level_neg_mesh_coh_bridge_perc', 'mt_coh_bridge_perc_mesh_w_level_neg',...
    'mt_bl_w_by_interp_coh_interp_grid_wneg', ...
    'mt_w_by_interp_coh_interp_grid_wneg', 'mt_w_by_interp_coh_interp_grid_wpos', 'mt_coh_w_perc_ratio_wneg'});
[ar_interp_c_grid, ar_interp_coh_grid, ar_a_meshk, ar_k_mesha, ...
    mt_interp_coh_grid_mesh_z, mt_z_mesh_coh_interp_grid, ...
    mt_interp_coh_grid_mesh_w_perc,...
    mt_w_level_neg_mesh_coh_bridge_perc, mt_coh_bridge_perc_mesh_w_level_neg, ...
    mt_bl_w_by_interp_coh_interp_grid_wneg, ...
    mt_w_by_interp_coh_interp_grid_wneg, mt_w_by_interp_coh_interp_grid_wpos, mt_coh_w_perc_ratio_wneg] ...
        = params_group{:};

params_group = values(armt_map, {'mt_coh_wkb', 'mt_z_mesh_coh_wkb'});
[mt_coh_wkb, mt_z_mesh_coh_wkb] = params_group{:};

% armt_map
% Formal choice Menu/Grid and Interest Rate Menu/Grid
params_group = values(armt_map, {'ar_forbrblk_r', 'ar_forbrblk'});
[ar_forbrblk_r, ar_forbrblk] = params_group{:};

% func_map
params_group = values(func_map, {'f_util_log', 'f_util_crra', 'f_cons'});
[f_util_log, f_util_crra, f_cons] = params_group{:};

% param_map
params_group = values(param_map, {'it_z_n', 'fl_crra', 'fl_beta', ...
    'fl_nan_replace', 'fl_c_min', 'bl_bridge', 'bl_default', 'fl_default_wprime'});
[it_z_n, fl_crra, fl_beta, fl_nan_replace, fl_c_min, bl_bridge, bl_default, fl_default_wprime] = params_group{:};
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

% collect optimal borrowing formal and informal choices
% mt_pol_b_with_r: cost to t+1 consumption from borrowing in t
mt_pol_b_with_r = zeros(length(ar_interp_coh_grid),length(ar_z));
mt_pol_b_bridge = zeros(length(ar_interp_coh_grid),length(ar_z));
mt_pol_inf_borr_nobridge = zeros(length(ar_interp_coh_grid),length(ar_z));
mt_pol_for_borr = zeros(length(ar_interp_coh_grid),length(ar_z));
mt_pol_for_save = zeros(length(ar_interp_coh_grid),length(ar_z));

% We did not need these in ff_oz_vf or ff_oz_vf_vec
% see
% <https://fanwangecon.github.io/M4Econ/support/speed/partupdate/fs_u_c_partrepeat_main.html
% fs_u_c_partrepeat_main> for why store using cells.
cl_u_c_store = cell([it_z_n, 1]);
cl_c_valid_idx = cell([it_z_n, 1]);
cl_w_kstar_interp_z = cell([it_z_n, 1]);
for it_z_i = 1:length(ar_z)
    cl_w_kstar_interp_z{it_z_i} = zeros([length(ar_w_perc), length(ar_interp_coh_grid)]) - 1;
end

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
    fl_u_cmin = f_util_log(fl_c_min);
else
    ar_interp_u_of_c_grid = f_util_crra(ar_interp_c_grid);
    fl_u_cmin = f_util_crra(fl_c_min);
end
ar_interp_u_of_c_grid(ar_interp_c_grid <= fl_c_min) = fl_u_cmin;

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
    % This is the same as <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbz/solve/html/ff_ipwkbz_vf_vecsv.html
    % ff_ipwkbz_vf_vecsv>. For the FIBS problem, the cash-on-hand
    % interpolation grid stays the same, and the shock grid stays the same
    % as well. The results will not be the same, for example, the coh_grid
    % max is the max of reachable cash-on-hand levels (min is however just
    % the borrowing bound).

    % Generate Interpolant for v(coh,z)
    f_grid_interpolant_value = griddedInterpolant(...
        mt_z_mesh_coh_interp_grid', mt_interp_coh_grid_mesh_z', mt_val_cur', 'linear', 'nearest');

    % Interpolate for v(coh(k(w,z),b(w,z),z),z)
    mt_val_wkb_interpolated = f_grid_interpolant_value(mt_z_mesh_coh_wkb, mt_coh_wkb);

    %% Solve Second Stage Problem k*(w,z)
    % This is again the same as <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbz/solve/html/ff_ipwkbz_vf_vecsv.html
    % ff_ipwkbz_vf_vecsv>. But the output matrix sizes are different.
    % Previously, they were (length(ar_w_level)) by (length(ar_z)). Now
    % have this thing which is stored (length(ar_w_level_full)) by
    % (length(ar_z)). _ar_w_level_full_ includes not just different levels
    % of _ar_w_level_, but also repeats the elements of _ar_w_level_ that
    % are < 0 by _it_coh_bridge_perc_n_ times, starting with what
    % corresponds to 100 percent of w should go to cover bridge loan, until
    % 0 percent for w < 0, which then proceeds to w > 0. So the last
    % segment of _ar_w_level_full_ is the same as ar_w_level:
    % ar_w_level_full((end-length(ar_w_level)+1):end) = ar_w_level.

    support_map('bl_graph_evf') = false;
    if (it_iter == (it_maxiter_val + 1))
        support_map('bl_graph_evf') = bl_graph_evf;
    end
    bl_input_override = true;
    [mt_ev_condi_z_max, ~, mt_ev_condi_z_max_kp, ~] = ...
        ff_ipwkbz_fibs_evf(mt_val_wkb_interpolated, param_map, support_map, armt_map, bl_input_override);

    %% Solve First Stage Problem w*(z) given k*(w,z)
    % Refer to
    % <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkbz/solve/html/ff_ipwkbz_fibs_vf_vecsv.html
    % ff_ipwkbz_fibs_vf_vecsv> where the problem was solved without formal and
    % informal choices that allow for bridge loans to see line by line how
    % code differ. Some of the comments from that file are not here to save
    % space. Comments here address differences and are specific to formal
    % and informal choices.

    % loop 1: over exogenous states
    for it_z_i = 1:length(ar_z)

        %% A. Interpolate FULL to get k*(coh_level, w_perc, z), b*(k,w) based on k*(coh_perc, w_level)
        % additionally, Interpolate FULL EV(k*(coh_level, w_perc, z), w -
        % b*|z) based on EV(k*(coh_perc, w_level))
        %
        % we solved the second period problem in ff_ipwkbz_fibs_fibs_evf.m
        % above. To use results, we need to interpolate in the following
        % way to obtain *mt_w_kstar_interp_z* as well as
        % *mt_ev_condi_z_max_interp_z*:
        %
        % # Interp STG1A: for $w > 0$, 1D interpolate over w level, given z
        % # Interp STG1B: for $w < 0$, 2D interpolate over w level and coh
        % perceng, given z
        %

        % 1. Negative Elements of w grid expanded by w percentages for bridge
        ar_bl_w_level_full_neg = (ar_w_level_full < 0);

        % 2. Current Positve w and negative w optimal k choices
        it_wneg_mt_row = sum(ar_bl_w_level_full_neg)/length(ar_coh_bridge_perc);
        % for mt_ev_condi_z_max_kp
        ar_ev_condi_z_max_kp_wpos = mt_ev_condi_z_max_kp(~ar_bl_w_level_full_neg, it_z_i)';
        ar_ev_condi_z_max_kp_wneg = mt_ev_condi_z_max_kp(ar_bl_w_level_full_neg, it_z_i)';
        mt_ev_condi_z_max_kp_wneg = reshape(ar_ev_condi_z_max_kp_wneg, [it_wneg_mt_row, length(ar_coh_bridge_perc)]);
        % for mt_ev_condi_z_max
        ar_ev_condi_z_max_wpos = mt_ev_condi_z_max(~ar_bl_w_level_full_neg, it_z_i)';
        ar_ev_condi_z_max_wneg = mt_ev_condi_z_max(ar_bl_w_level_full_neg, it_z_i)';
        mt_ev_condi_z_max_wneg = reshape(ar_ev_condi_z_max_wneg, [it_wneg_mt_row, length(ar_coh_bridge_perc)]);

        % 2. Interp STG1A for w > 0
        ar_w_level_full_pos = ar_w_level_full(~ar_bl_w_level_full_neg);
        % Interpolation for mt_ev_condi_z_max_kp
        f_interpolante_w_level_pos_kstar_z = griddedInterpolant(ar_w_level_full_pos, ar_ev_condi_z_max_kp_wpos, 'linear', 'nearest');
        mt_w_kstar_interp_z_wpos = f_interpolante_w_level_pos_kstar_z(mt_w_by_interp_coh_interp_grid_wpos(:));
        mt_w_astar_interp_z_wpos = mt_w_by_interp_coh_interp_grid_wpos(:) - mt_w_kstar_interp_z_wpos;
        % Interpolation for mt_ev_condi_z_max
        f_interpolante_w_level_pos_ev_z = griddedInterpolant(ar_w_level_full_pos, ar_ev_condi_z_max_wpos, 'linear', 'nearest');
        mt_w_ev_interp_z_wpos = f_interpolante_w_level_pos_ev_z(mt_w_by_interp_coh_interp_grid_wpos(:));

        % 3. Interp STG1B for w <= 0
        if (bl_bridge)
            % Interpolation for mt_ev_condi_z_max_kp
            f_interpolante_w_level_neg_kstar_z = griddedInterpolant(...
                mt_coh_bridge_perc_mesh_w_level_neg', mt_w_level_neg_mesh_coh_bridge_perc', ...
                mt_ev_condi_z_max_kp_wneg', 'linear', 'nearest');
            mt_w_kstar_interp_z_wneg = f_interpolante_w_level_neg_kstar_z(mt_coh_w_perc_ratio_wneg(:), mt_w_by_interp_coh_interp_grid_wneg(:));
            mt_w_astar_interp_z_wneg = mt_w_by_interp_coh_interp_grid_wneg(:) - mt_w_kstar_interp_z_wneg;
            % Interpolation for mt_ev_condi_z_max
            f_interpolante_w_level_neg_ev_z = griddedInterpolant(...
                mt_coh_bridge_perc_mesh_w_level_neg', mt_w_level_neg_mesh_coh_bridge_perc', ...
                mt_ev_condi_z_max_wneg', 'linear', 'nearest');
            mt_w_ev_interp_z_wneg = f_interpolante_w_level_neg_ev_z(mt_coh_w_perc_ratio_wneg(:), mt_w_by_interp_coh_interp_grid_wneg(:));
        else
            ar_w_level_full_neg = ar_w_level_full(ar_bl_w_level_full_neg);
            % Interpolation for mt_ev_condi_z_max_kp
            f_interpolante_w_level_neg_kstar_z = griddedInterpolant(ar_w_level_full_neg, ar_ev_condi_z_max_kp_wneg, 'linear', 'nearest');
            mt_w_kstar_interp_z_wneg = f_interpolante_w_level_neg_kstar_z(mt_w_by_interp_coh_interp_grid_wneg(:));
            mt_w_astar_interp_z_wneg = mt_w_by_interp_coh_interp_grid_wneg(:) - mt_w_kstar_interp_z_wneg;
            % Interpolation for mt_ev_condi_z_max
            f_interpolante_w_level_neg_ev_z = griddedInterpolant(ar_w_level_full_neg, ar_ev_condi_z_max_wneg, 'linear', 'nearest');
            mt_w_ev_interp_z_wneg = f_interpolante_w_level_neg_ev_z(mt_w_by_interp_coh_interp_grid_wneg(:));
        end

        % 4. Combine positive and negative aggregate savings matrix
        % check: mt_w_by_interp_coh_interp_grid vs mt_w_astar_interp_z + mt_w_kstar_interp_z
        % combine for mt_ev_condi_z_max_kp
        mt_w_kstar_interp_z = zeros(size(mt_bl_w_by_interp_coh_interp_grid_wneg));
        mt_w_kstar_interp_z(~mt_bl_w_by_interp_coh_interp_grid_wneg) = mt_w_kstar_interp_z_wpos;
        mt_w_kstar_interp_z(mt_bl_w_by_interp_coh_interp_grid_wneg)  = mt_w_kstar_interp_z_wneg;
        mt_w_astar_interp_z = zeros(size(mt_bl_w_by_interp_coh_interp_grid_wneg));
        mt_w_astar_interp_z(~mt_bl_w_by_interp_coh_interp_grid_wneg) = mt_w_astar_interp_z_wpos;
        mt_w_astar_interp_z(mt_bl_w_by_interp_coh_interp_grid_wneg)  = mt_w_astar_interp_z_wneg;
        % combine for mt_ev_condi_z_max
        mt_ev_condi_z_max_interp_z = zeros(size(mt_bl_w_by_interp_coh_interp_grid_wneg));
        mt_ev_condi_z_max_interp_z(~mt_bl_w_by_interp_coh_interp_grid_wneg) = mt_w_ev_interp_z_wpos;
        mt_ev_condi_z_max_interp_z(mt_bl_w_by_interp_coh_interp_grid_wneg)  = mt_w_ev_interp_z_wneg;

        % 5. changes in w_perc kstar choices
        mt_w_kstar_diff_idx = (cl_w_kstar_interp_z{it_z_i} ~= mt_w_kstar_interp_z);

        %% B. Calculate UPDATE u(c) Update: u(c(coh_level, w_perc)) given k*_interp, b*_interp
        ar_c = f_cons(mt_interp_coh_grid_mesh_w_perc(mt_w_kstar_diff_idx), ...
                      mt_w_astar_interp_z(mt_w_kstar_diff_idx), ...
                      mt_w_kstar_interp_z(mt_w_kstar_diff_idx));

        ar_it_c_valid_idx = (ar_c <= fl_c_min);
        % EVAL current utility: N by N, f_util defined earlier
        ar_utility_update = f_grid_interpolant_spln(ar_c);

        % Update Storage
        if (it_iter == 1)
            cl_u_c_store{it_z_i} = reshape(ar_utility_update, [length(ar_w_perc), length(ar_interp_coh_grid)]);
            cl_c_valid_idx{it_z_i} = reshape(ar_it_c_valid_idx, [length(ar_w_perc), length(ar_interp_coh_grid)]);
        else
            cl_u_c_store{it_z_i}(mt_w_kstar_diff_idx) = ar_utility_update;
            cl_c_valid_idx{it_z_i}(mt_w_kstar_diff_idx) = ar_it_c_valid_idx;
        end
        cl_w_kstar_interp_z{it_z_i} = mt_w_kstar_interp_z;

        %% D. Compute FULL U(coh_level, w_perc, z) over all w_perc
        mt_utility = cl_u_c_store{it_z_i} + fl_beta*mt_ev_condi_z_max_interp_z;

        % Index update
        % using the method below is much faster than index replace
        % see <https://fanwangecon.github.io/M4Econ/support/speed/index/fs_subscript.html fs_subscript>
        mt_it_c_valid_idx = cl_c_valid_idx{it_z_i};
        % Default or Not Utility Handling
        if (bl_default)
            % if default: only today u(cmin), transition out next period, debt wiped out
            fl_v_default = fl_u_cmin + fl_beta*f_interpolante_w_level_pos_ev_z(fl_default_wprime);
            mt_utility = mt_utility.*(~mt_it_c_valid_idx) + fl_v_default*(mt_it_c_valid_idx);
        else
            % if default is not allowed: v = u(cmin)
            mt_utility = mt_utility.*(~mt_it_c_valid_idx) + fl_nan_replace*(mt_it_c_valid_idx);
        end

        % percentage algorithm does not have invalid (check to make sure
        % min percent is not 0 in ffs_ipwkbz_fibs_get_funcgrid.m)
        % mt_utility = mt_utility.*(~mt_it_c_valid_idx) + fl_u_neg_c*(mt_it_c_valid_idx);

        %% E. Optimize Over Choices: max_{w_perc} U(coh_level, w_perc, z)
        % Optimization: remember matlab is column major, rows must be
        % choices, columns must be states
        % <https://en.wikipedia.org/wiki/Row-_and_column-major_order COLUMN-MAJOR>
        [ar_opti_val_z, ar_opti_idx_z] = max(mt_utility);

        % Generate Linear Opti Index
        [it_choies_n, it_states_n] = size(mt_utility);
        ar_add_grid = linspace(0, it_choies_n*(it_states_n-1), it_states_n);
        ar_opti_linear_idx_z = ar_opti_idx_z + ar_add_grid;

        ar_opti_aprime_z = mt_w_astar_interp_z(ar_opti_linear_idx_z);
        ar_opti_kprime_z = mt_w_kstar_interp_z(ar_opti_linear_idx_z);
        ar_opti_c_z = f_cons(ar_interp_coh_grid, ar_opti_aprime_z, ar_opti_kprime_z);

        % Handle Default is optimal or not
        if (bl_default)
            % if defaulting is optimal choice, at these states, not required
            % to default, non-default possible, but default could be optimal
            fl_default_opti_kprime = f_interpolante_w_level_pos_kstar_z(fl_default_wprime);
            ar_opti_aprime_z(ar_opti_c_z <= fl_c_min) = fl_default_wprime - fl_default_opti_kprime;
            ar_opti_kprime_z(ar_opti_c_z <= fl_c_min) = fl_default_opti_kprime;
        else
            % if default is not allowed, then next period same state as now
            % this is absorbing state, this is the limiting case, single
            % state space point, lowest a and lowest shock has this.
            ar_opti_aprime_z(ar_opti_c_z <= fl_c_min) = min(ar_a_meshk);
            ar_opti_kprime_z(ar_opti_c_z <= fl_c_min) = min(ar_k_mesha);
        end

        %% F. Store Results
        mt_val(:,it_z_i) = ar_opti_val_z;
        mt_pol_a(:,it_z_i) = ar_opti_aprime_z;
        mt_pol_k(:,it_z_i) = ar_opti_kprime_z;
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

%% Process Optimal Choices 1: Formal and Informal Choices

result_map = containers.Map('KeyType','char', 'ValueType','any');
result_map('mt_val') = mt_val;
result_map('mt_pol_idx') = mt_pol_idx;

% Find optimal Formal Informal Choices. Could have saved earlier, but was
% wasteful of resources
for it_z_i = 1:length(ar_z)
    for it_coh_interp_j = 1:length(ar_interp_coh_grid)

        fl_coh = mt_interp_coh_grid_mesh_z(it_coh_interp_j, it_z_i);
        fl_a_opti = mt_pol_a(it_coh_interp_j, it_z_i);

        % call formal and informal function.
        [fl_max_c, fl_opti_b_bridge, fl_opti_inf_borr_nobridge, fl_opti_for_borr, fl_opti_for_save] = ...
            ffs_fibs_min_c_cost_bridge(fl_a_opti, fl_coh, param_map, support_map, armt_map, func_map);

        % store savings and borrowing formal and inf optimal choices
        mt_pol_b_with_r(it_coh_interp_j,it_z_i) = fl_max_c;
        mt_pol_b_bridge(it_coh_interp_j,it_z_i) = fl_opti_b_bridge;
        mt_pol_inf_borr_nobridge(it_coh_interp_j,it_z_i) = fl_opti_inf_borr_nobridge;
        mt_pol_for_borr(it_coh_interp_j,it_z_i) = fl_opti_for_borr;
        mt_pol_for_save(it_coh_interp_j,it_z_i) = fl_opti_for_save;

    end
end

%% Process Optimal Choices 2: Store a, k, c, coh Results
%
% # *mt_interp_coh_grid_mesh_z*: Cash-on-hand period _t_.
% # *mt_pol_a*: Safe asset choice, principles only for ipwkbz_fibs
% # *cl_mt_pol_a_principleonly*: mt_pol_a is stored in
% cl_mt_pol_a_principle only. This is a shortcut because we need to keep
% cl_mt_pol_a for mt_pol_b_with_r for the _ds_ code.
% # *cl_mt_pol_a*: stores _mt_pol_b_with_r_ which has principles and
% interest rates, to be used with _ds_ code.
% # *mt_pol_a*: Safe asset choice, principles only for ipwkbz_fibs
% # *cl_mt_pol_k*: Risky asset choice, principles only for ipwkbz_fibs
% # *cl_mt_pol_c*: Consumption in _t_ given choices.
% # *cl_pol_b_with_r*: Consumption cost of _mt_pol_a_ in _t+1_, given the
% formal and informal choices that are optimal to minimize this consumption
% cost.
%

result_map('cl_mt_val') = {mt_val, zeros(1)};
result_map('cl_mt_coh') = {mt_interp_coh_grid_mesh_z, zeros(1)};
result_map('cl_mt_pol_k') = {mt_pol_k, zeros(1)};
result_map('cl_mt_pol_c') = {f_cons(mt_interp_coh_grid_mesh_z, mt_pol_a, mt_pol_k), zeros(1)};

result_map('cl_mt_pol_a') = {mt_pol_b_with_r, zeros(1)};
result_map('cl_mt_pol_a_principleonly') = {mt_pol_a, zeros(1)};

%% Process Optimal Choices 3: Store Formal and Informal Choices
result_map('cl_mt_pol_b_bridge') = {mt_pol_b_bridge, zeros(1)};
result_map('cl_mt_pol_inf_borr_nobridge') = {mt_pol_inf_borr_nobridge, zeros(1)};
result_map('cl_mt_pol_for_borr') = {mt_pol_for_borr, zeros(1)};
result_map('cl_mt_pol_for_save') = {mt_pol_for_save, zeros(1)};

%% Process Optimal Choices 4: List of Variable Names to be processed by distributional codes
% this list is needed for the ds codes to generate distribution,
% distributional statistcs will be computed for elements in the list here.

result_map('ar_st_pol_names') = ...
    ["cl_mt_val", "cl_mt_coh", "cl_mt_pol_a", "cl_mt_pol_k", "cl_mt_pol_c", "cl_mt_pol_a_principleonly", ...
    "cl_mt_pol_b_bridge", "cl_mt_pol_inf_borr_nobridge", "cl_mt_pol_for_borr", "cl_mt_pol_for_save"];

% Get Discrete Choice Outcomes
result_map = ffs_fibs_identify_discrete(result_map, bl_input_override);

%% End Timer and Profile

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

%% Post Solution Graph and Table Generation

if (bl_post)
    bl_input_override = true;
    result_map('ar_val_diff_norm') = ar_val_diff_norm(1:it_iter_last);
    result_map('ar_pol_diff_norm') = ar_pol_diff_norm(1:it_iter_last);
    result_map('mt_pol_perc_change') = mt_pol_perc_change(1:it_iter_last, :);

    armt_map('mt_coh_wkb_ori') = mt_coh_wkb;
    armt_map('ar_a_meshk_ori') = ar_a_meshk;
    armt_map('ar_k_mesha_ori') = ar_k_mesha;
    
    armt_map('mt_coh_wkb') = mt_interp_coh_grid_mesh_z;
    armt_map('it_ameshk_n') = length(ar_interp_coh_grid);
    armt_map('ar_a_meshk') = mt_interp_coh_grid_mesh_z(:,1);
    armt_map('ar_k_mesha') = zeros(size(mt_interp_coh_grid_mesh_z(:,1)) + 0);

    % Standard AZ graphs
    result_map = ff_akz_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

    % Graphs for results_map with FIBS contents
    armt_map('ar_a') = ar_interp_coh_grid;
    result_map = ff_az_fibs_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);

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
