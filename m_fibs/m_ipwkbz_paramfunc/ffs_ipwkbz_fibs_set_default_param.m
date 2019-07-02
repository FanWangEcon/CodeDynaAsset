%% Set Model Parameters (Interpolated + Percentage + Risky + Safe Asset + FIBS)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [param_map, support_map] = ffs_ipwkbz_fibs_set_default_param(varargin)
%% FFS_IPWKBZ_FIBS_SET_DEFAULT_PARAM setting model default parameters

%% Default

it_subset = 0;
bl_display_defparam = false;
default_params = {it_subset bl_display_defparam};
[default_params{1:length(varargin)}] = varargin{:};
[it_subset, bl_display_defparam] = default_params{:};

%% Setting param_map container

param_map = containers.Map('KeyType','char', 'ValueType','any');

% Preferences
param_map('fl_crra') = 1.5;
param_map('fl_beta') = 0.94;

% Production Function
% Productivity Shock Parameters
param_map('it_z_n') = 15;
param_map('fl_z_mu') = 0;
param_map('fl_z_rho') = 0.8;
param_map('fl_z_sig') = 0.2;

% Borrowing
% fl_default_aprime is the next period asset level
% households face if they default.
param_map('fl_b_bd') = -20; % borrow bound, = 0 if save only
param_map('fl_default_wprime') = 0; % wprime not a prime

% Borrowing Setting 1: Default Allowed, Bridge True, bl_rollover does not matter
% Borrowing Setting 2: Default Allowed, Bridge False, bl_rollover matter
param_map('bl_default') = true; % if borrowing is default allowed
param_map('bl_bridge') = false;
param_map('bl_rollover') = true;

% CD Production Function Parameters
param_map('fl_Amean') = 1;
param_map('fl_alpha') = 0.36;
param_map('fl_delta') = 0.08;

% Prices
% shock is on k, not on labor, fl_w is fixed wage income
param_map('fl_w') = 1.28*0.3466; % min(z*w) from benchmark az model

% formal informal parameters
% fl_for_br_block are the formal borrowing grid block sizes.
param_map('fl_r_fsv') = 0.025;
param_map('fl_r_inf') = 0.045;
param_map('fl_r_inf_bridge') = 0.045;
param_map('fl_r_fbr') = 0.035;
param_map('bl_b_is_principle') = true;
% see: ffs_for_br_block.m
param_map('st_forbrblk_type') = 'seg3';
param_map('fl_forbrblk_brmost') = -19;
param_map('fl_forbrblk_brleast') = -1;
param_map('fl_forbrblk_gap') = -1.5;

% Minimum Consumption, utility lower bound. The cmin parameter and
% fl_nan_replace parameter have no effects on value function, just for
% resetting invalid choice grid values. fl_nan_replace reset invalid k
% choice given w. fl_c_min resets invalid consumption levels due to w
% choices that are invalid. But this is the case when fl_w > 0.
param_map('fl_c_min') = 0.001;
param_map('fl_nan_replace') = -9999;

% Asset Grids
% Toal savings aggregate grid (see discussion on top). 35 points picked for
% for this problem w_max is overall for everyone, but each individual coh
% levle has associated w_max. Unlike before, when we had it_w_n, now we
% have it_w_perc_n which is how many percentage grid points to have. We
% also now include fl_w_interp_grid_gap, the grip gap for interpolation.
param_map('fl_w_min') = param_map('fl_b_bd'); % but b_bd overrides this
param_map('fl_w_max') = 50;
param_map('it_w_perc_n') = 50;

% Risky Capital Asset Vector
% see graph below for how it looks graphically
% in principle keep it_k_n the same as it_w_n to have equi-distance points
% in triangle choice grid. For the ipwkz this k_max is the overall max, not
% the individual max given coh. we have it_ak_perc_n rather than it_ak_n
% for percentage.
param_map('fl_k_min') = 0;
param_map('fl_k_max') = (param_map('fl_w_max') - param_map('fl_b_bd'));
param_map('it_ak_perc_n') = param_map('it_w_perc_n'); % grid for a and k the same

% Percentage of w that is not for bridge loan, when param_map('bl_bridge') = false
% ar_coh_bridge_perc = [1]
param_map('it_coh_bridge_perc_n') = param_map('it_w_perc_n');

% Interpolation
% fl_coh_interp_grid_gap controls the number of coh points at which to solve the model
% it_c_interp_grid_gap determines the gap between consumption terpolation
% points. For consumption interpolation 10^-4 is extremely accurate, there
% should be no perceptible differences in value and policy functions when the
% it_c_interp_grid_gap <= 0.001 compared to actual evaluation. Also include
% now fl_w_interp_grid_gap above, which is for interpolation over w.
param_map('fl_coh_interp_grid_gap') = 0.1;
% param_map('it_coh_interp_n') = 500;
param_map('it_c_interp_grid_gap') = 10^-4;
% Interpolation gap second stage w
% previously only it_w_n, now two grids control w it_w_perc_n for 2nd stage
% fl_w_interp_grid_gap for first stage, make them the same length for
% default.
param_map('fl_w_interp_grid_gap') = 0.1;
% param_map('fl_w_interp_grid_gap') = (param_map('fl_w_max') - param_map('fl_w_min'))/param_map('it_w_perc_n');

% Solution Accuracy
param_map('it_maxiter_val') = 250;
param_map('it_maxiter_dist') = 250;
param_map('it_trans_power_dist') = 250;
param_map('st_analytical_stationary_type') = 'eigenvector'; % could be eigenvector, projection, power

param_map('fl_tol_val') = 10^-5;
param_map('fl_tol_pol') = 10^-5;
param_map('fl_tol_dist') = 10^-5;
param_map('it_tol_pol_nochange') = 25; % number of iterations where policy does not change

%% Setting support_map container

support_map = containers.Map('KeyType','char', 'ValueType','any');

% root directory
[st_root_path] = preamble(false);
st_matimg_path_root = [st_root_path '/m_fibs/'];
support_map('st_matimg_path_root') = st_matimg_path_root;

% timer
support_map('bl_time') = true;

% Print Controls
support_map('bl_display') = true;
support_map('bl_display_dist') = false;
support_map('it_display_every') = 5; % how often to print results

% Profile Controls
support_map('bl_profile') = false;
support_map('bl_profile_dist') = false; % distribution profile
support_map('st_profile_path') = [st_matimg_path_root '/m_ipwkbz_solve/profile/'];
support_map('st_profile_prefix') = [''];
support_map('st_profile_name_main') = ['_default'];
support_map('st_profile_suffix') = ['_p' num2str(it_subset)];

support_map('bl_post') = false;

% Final Print
support_map('bl_display_final') = false; % print finalized results
support_map('bl_display_final_dist') = false; % print finalized results
support_map('it_display_final_rowmax') = 100; % max row to print (states/iters)
support_map('it_display_final_colmax') = 12; % max col to print (shocks)

% Mat File Controls
support_map('bl_mat') = false;
support_map('st_mat_path') = [st_matimg_path_root '/m_ipwkbz_solve/mat/'];
support_map('st_mat_prefix') = [''];
support_map('st_mat_name_main') = ['_default'];
support_map('st_mat_suffix') = ['_p' num2str(it_subset)];

% Graphing Controls
support_map('bl_graph') = false;
support_map('bl_graph_onebyones') = true;
support_map('bl_graph_val') = true;
support_map('bl_graph_pol_lvl') = true;
support_map('bl_graph_pol_pct') = true;
support_map('bl_graph_coh_t_coh') = true;
support_map('bl_graph_discrete') = true;
% Formal Informal Specific Graphs
support_map('bl_graph_forinf_discrete') = true;
support_map('bl_graph_forinf_pol_lvl') = true;
support_map('bl_graph_forinf_pol_pct') = true;


% Image Saving Controls (given graphing)
support_map('st_title_prefix') = '';
support_map('bl_img_save') = false;
support_map('st_img_path') = [st_matimg_path_root '/m_ipwkbz_solve/img/'];
support_map('st_img_prefix') = [''];
support_map('st_img_name_main') = ['_default'];
support_map('st_img_suffix') = ['_p' num2str(it_subset) '.png'];

% Sub-function graphing controls
support_map('bl_graph_funcgrids') = false;
support_map('bl_graph_funcgrids_detail') = false;
support_map('bl_display_funcgrids') = false;
support_map('bl_graph_evf') = false;
support_map('bl_display_evf') = false;
support_map('bl_display_minccost') = false;
support_map('bl_display_infbridge') = false;


%% Subset Options
%
% # it_subset = 1 is basic invoke quick test
% # it_subset = 2 is main invoke
% # it_subset = 3 is profiling invoke
% # it_subset = 4 is matlab publish.
%

if (ismember(it_subset, [1,2,3,4]))
    if (ismember(it_subset, [1]))
        % TEST quick
        param_map('it_w_perc_n') = 20;
        param_map('it_ak_perc_n') = param_map('it_w_perc_n');
        param_map('it_z_n') = 3;

        param_map('fl_coh_interp_grid_gap') = 0.25;
        param_map('it_c_interp_grid_gap') = 0.001;
        param_map('fl_w_interp_grid_gap') = 1;

        param_map('it_maxiter_val') = 50;
        param_map('it_tol_pol_nochange') = 1000;
        support_map('bl_display') = true;
        support_map('it_display_every') = 1;
    end
    if (ismember(it_subset, [2, 4]))
        % close figures
        close all;
        % Main Run
        support_map('bl_time') = true;
        support_map('bl_display') = true;
        support_map('it_display_every') = 5;

        support_map('bl_post') = true;
        support_map('bl_display_final') = true;
        support_map('bl_mat') = false;
        support_map('bl_graph') = true;
        support_map('bl_graph_onebyones') = false;
        support_map('bl_img_save') = true;
        if (ismember(it_subset, [4]))
            support_map('bl_time') = false;
            support_map('bl_display') = false;
            support_map('bl_graph_onebyones') = true;
            support_map('bl_img_save') = false;
        end
    end
    if (ismember(it_subset, [3]))
        % Profile run
        support_map('bl_profile') = true;
        support_map('bl_display') = false; % don't print
        support_map('bl_time') = true;
    end
end


%% Subset Options for Distribution solutions
%
% # it_subset = 5 is basic invoke quick test
% # it_subset = 6 is invoke full test
% # it_subset = 7 is profiling invoke
% # it_subset = 8 is matlab publish
% # it_subset = 9 is invoke operational (only final stats) and coh graph
%

if (ismember(it_subset, [5,6,7,8,9]))
    if (ismember(it_subset, [5]))
        % TEST quick (need to enough to have distribution)

        param_map('it_w_perc_n') = 100;
        param_map('it_ak_perc_n') = param_map('it_w_perc_n');
        param_map('it_z_n') = 7;

        param_map('fl_coh_interp_grid_gap') = 0.25;
        param_map('it_c_interp_grid_gap') = 0.001;
        param_map('fl_w_interp_grid_gap') = 1;

        param_map('it_maxiter_val') = 50;
        param_map('it_tol_pol_nochange') = 1000;
        support_map('bl_display') = true;
        support_map('it_display_every') = 1;

        support_map('bl_display_dist') = true;
    end
    if (ismember(it_subset, [6, 8, 9]))
        % close all
        close all;
        % Main Run
        support_map('bl_time') = true;
        support_map('bl_display') = false;
        support_map('bl_display_dist') = true;
        support_map('it_display_every') = 20;

        support_map('bl_post') = true;
        support_map('bl_display_final_dist') = true;
        support_map('bl_mat') = false;
        support_map('bl_graph') = true;
        support_map('bl_graph_onebyones') = false;
        support_map('bl_img_save') = true;

        % do not generate all graphs when solving for distribution
        support_map('bl_graph_val') = false;
        support_map('bl_graph_pol_lvl') = false;
        support_map('bl_graph_pol_pct') = false;
        support_map('bl_graph_coh_t_coh') = true;
        support_map('bl_graph_forinf_discrete') = false;
        support_map('bl_graph_forinf_pol_lvl') = false;
        support_map('bl_graph_forinf_pol_pct') = true;

        if (ismember(it_subset, [8, 9]))
            support_map('bl_time') = false;
            support_map('bl_display') = false;
            support_map('bl_display_dist') = false;
            support_map('bl_graph_onebyones') = true;
            support_map('bl_img_save') = false;
            if (ismember(it_subset, [9]))
                support_map('bl_graph_coh_t_coh') = false;
                support_map('bl_graph_forinf_pol_pct') = false;
            end
        end

    end
    if (ismember(it_subset, [7]))
        % Profile run
        support_map('bl_profile_dist') = true;
        support_map('bl_display') = false; % don't print
        support_map('bl_display_dist') = false; % don't print
        support_map('bl_time') = true;
    end
end


%% Display

if (bl_display_defparam)
    disp('param_map');
    disp(param_map);
    param_map_keys = keys(param_map);
    param_map_vals = values(param_map);
    for i = 1:length(param_map)
        st_display = strjoin(['pos =' num2str(i) '; key =' string(param_map_keys{i}) '; val =' string(param_map_vals{i})]);
        disp(st_display);
    end

    disp('support_map')
    disp(support_map);
    param_map_keys = keys(support_map);
    param_map_vals = values(support_map);
    for i = 1:length(support_map)
        st_display = strjoin(['pos =' num2str(i) '; key =' string(param_map_keys{i}) '; val =' string(param_map_vals{i})]);
        disp(st_display);
    end
end

end
