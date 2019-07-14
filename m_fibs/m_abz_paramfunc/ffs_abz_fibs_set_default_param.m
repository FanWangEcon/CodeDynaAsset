%% Set Model Parameters (ABZ FIBS)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [param_map, support_map] = ffs_abz_fibs_set_default_param(varargin)
%% FFS_ABZ_FIBS_SET_DEFAULT_PARAM setting model default parameters
% two groups of default parameters stored in container maps
%
% @param it_subset integer default parameter control subsetting. it_subset = 1 is
% basic invoke quick test. it_subset = 2 is main invoke. it_subset = 3 is
% profiling invoke. it_subset = 4 is matlab publish.
%
% @param bl_display_defparam boolean local printing
%
% @return param_map container parameters needed for solving the model
%
% @return support_map container programming control parameters like to graph to print etc
%
% @example
%
%   it_param_set = 1;
%   [param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);
%

%% Default

it_subset = 4;
if (isempty(varargin))
    bl_display_defparam = true;
else
    bl_display_defparam = false;
end
default_params = {it_subset bl_display_defparam};
[default_params{1:length(varargin)}] = varargin{:};
[it_subset, bl_display_defparam] = default_params{:};

%% Setting param_map container

param_map = containers.Map('KeyType','char', 'ValueType','any');

param_map('st_model') = 'abz_fibs';

% Preferences
param_map('fl_crra') = 1.5;
param_map('fl_beta') = 0.94;

% Shock Parameters
param_map('it_z_n') = 15;
param_map('fl_z_mu') = 0;
param_map('fl_z_rho') = 0.8;
param_map('fl_z_sig') = 0.2;

% Borrowing
% fl_default_aprime is the next period asset level
% households face if they default.
param_map('fl_b_bd') = -20; % borrow bound, = 0 if save only
param_map('fl_default_aprime') = 0;

% Borrowing Setting 1: Default Allowed, Bridge True, bl_rollover does not matter
% Borrowing Setting 2: Default Allowed, Bridge False, bl_rollover matter
param_map('bl_default') = true; % if borrowing is default allowed
param_map('bl_bridge') = true;
param_map('bl_rollover') = true;

% Savings
param_map('fl_a_min') = 0; % if there is minimum savings requirement
param_map('fl_a_max') = 50;
param_map('bl_loglin') = false; % log lin threshold structure
param_map('fl_loglin_threshold') = 1; % dense points before 1
param_map('it_a_n') = 750;

% Prices
param_map('fl_w') = 1.28;

% formal informal parameters
% fl_for_br_block are the formal borrowing grid block sizes.
param_map('fl_r_fsv') = 0.025;
param_map('fl_r_inf') = 0.095;
param_map('fl_r_inf_bridge') = 0.095;
param_map('fl_r_fbr') = 0.065;
param_map('bl_b_is_principle') = false;
% see: ffs_for_br_block.m
param_map('st_forbrblk_type') = 'seg3';
param_map('fl_forbrblk_brmost') = -19;
param_map('fl_forbrblk_brleast') = -1;
param_map('fl_forbrblk_gap') = -1.5;

% Minimum Consumption, c_min is for default, when c < 0, replace utility
% with fl_nan_replace.
param_map('fl_c_min') = 0.02;
param_map('fl_nan_replace') = -99999;

% Solution Accuracy
param_map('it_maxiter_val') = 1000;
param_map('it_maxiter_dist') = 1000;
param_map('it_trans_power_dist') = 1000;
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
support_map('st_profile_path') = [st_matimg_path_root '/m_abz_solve/profile/'];
support_map('st_profile_prefix') = [''];
support_map('st_profile_name_main') = ['_default'];
support_map('st_profile_suffix') = ['_p' num2str(it_subset)];

support_map('bl_post') = false;
% Final Print
support_map('bl_display_final') = false; % print finalized results
support_map('bl_display_final_dist') = false; % print finalized results
support_map('bl_display_final_dist_detail') = false; % print finalized results
support_map('it_display_final_rowmax') = 100; % max row to print (states/iters)
support_map('it_display_final_colmax') = 12; % max col to print (shocks)
% Mat File Controls
support_map('bl_mat') = false;
support_map('st_mat_path') = [st_matimg_path_root '/m_abz_solve/mat/'];
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
support_map('st_img_path') = [st_matimg_path_root '/m_abz_solve/img/'];
support_map('st_img_prefix') = [''];
support_map('st_img_name_main') = ['_default'];
support_map('st_img_suffix') = ['_p' num2str(it_subset) '.png'];

% Sub-function graphing controls
support_map('bl_graph_funcgrids') = false;
support_map('bl_display_funcgrids') = false;
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
        param_map('it_a_n') = 25;
        param_map('it_z_n') = 3;
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
        param_map('it_a_n') = 100;
        param_map('it_z_n') = 7;
        param_map('it_maxiter_val') = 50;
        param_map('it_maxiter_dist') = 50;
        param_map('it_tol_pol_nochange') = 1000;
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
            support_map('bl_display_final_dist_detail') = true;
            support_map('bl_graph_onebyones') = true;
            support_map('bl_img_save') = false;
            if (ismember(it_subset, [9]))
                support_map('bl_display_final_dist_detail') = false;
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

%% Display

if (bl_display_defparam)
    fft_container_map_display(param_map);
    fft_container_map_display(support_map);
end

end
