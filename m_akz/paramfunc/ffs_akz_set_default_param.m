%% Set Model Parameters (Risky + Safe Asset)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [param_map, support_map] = ffs_akz_set_default_param(varargin)
%% FFS_AKZ_SET_DEFAULT_PARAM setting model default parameters
% two groups of default parameters stored in container maps. Explicitly
% solving for both a and k grids exponentially increase problem size, and
% we can not have state and choice grids that are as tight as for the az
% problem. Note that the defaults set it_w_n to have a far more limited
% number of points. For each w, we allow for combinations of a and k
% choices. Before, each aggregate savings choice only had one a associated
% with it.
%
% The biggest change from
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html
% ffs_az_set_default_param> is that now we have two asset choices. Previous
% we had the a grid for the cash-on-hand and z (az) model.
%
% Now a is w. There was only one asset before, so a was total savings. Now
% w is total savings w = a' + k'. Define choice grid in terms of the total
% savings, rather than a and k grids separately. If we define a and k with
% min and max separately, we will have a choice combo max(a) + max(k), this
% corresponds to one very high level of aggregate savings, but we are
% implicitly only allowing for dividing a and k evenly when we choose this.
% Implicitly, the higher the aggregate savings, the more limited the a and
% k combinations becomes when we define a and k grids separately.
% Geometrically, With aprime as the x-axis and kprime as the y-axis of
% choices, we need to define choices in terms of diagnoals (triangles)
% rather than vertical and horizontal lines (rectangles). See this
% graphically at the end of this page:
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html
% ffs_akz_get_funcgrid> which shows the triangular choice grids.
%
% With our triangular shape, with w_n total saving grid points, we will
% have : w_n*(w_n-1)/2 + w_n total grid of a and k jointly. Specifically,
% with the benchmark grid of w_n = 45 points, we have (45-1)*45/2 + 45 =
% 1035 total a and k joint grid points.
%
% @param it_subset integer default parameter control subsetting. it_subset
% = 1 is basic invoke quick test. it_subset = 2 is main invoke. it_subset =
% 3 is profiling invoke. it_subset = 4 is matlab publish.
%
% @param bl_display_defparam boolean local printing
%
% @return param_map container parameters needed for solving the model
%
% @return support_map container programming control parameters like to
% graph to print etc
%
% @example
%
%   it_param_set = 1; [param_map, support_map] =
%   ffs_akz_set_default_param(it_param_set);
%

%% Default

it_subset = 0;
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

% model name
param_map('st_model') = 'akz_wkz_iwkz';

% Preferences
param_map('fl_crra') = 1.5;
param_map('fl_beta') = 0.94;

% Production Function Productivity Shock Parameters
param_map('it_z_n') = 15;
param_map('fl_z_mu') = 0;
param_map('fl_z_rho') = 0.8;
param_map('fl_z_sig') = 0.2;

% CD Production Function Parameters
param_map('fl_Amean') = 1;
param_map('fl_alpha') = 0.36;
param_map('fl_delta') = 0.08;

% Prices shock is on k, not on labor, fl_w is fixed wage income
param_map('fl_w') = 1.28*0.3466; % min(z*w) from benchmark az model
param_map('fl_r_save') = 0.025;
param_map('fl_r_borr') = 0.025;

% Minimum Consumption, utility lower bound. The cmin parameter and
% fl_nan_replace parameter have no effects on value function, just for
% resetting invalid choice grid values. fl_nan_replace reset invalid k
% choice given w. fl_c_min resets invalid consumption levels due to w
% choices that are invalid. But this is the case when fl_w > 0.
param_map('fl_c_min') = 0.001;
param_map('fl_nan_replace') = -9999;

% Borrowing Related Parameters for AKBZ models the
% <https://github.com/FanWangEcon/CodeDynaAsset/tree/master/m_akz m_akz>
% models do not deal with default,
% <https://github.com/FanWangEcon/CodeDynaAsset/tree/master/m_akbz m_akbz>
% files do.
param_map('fl_b_bd') = 0; % borrow bound, = 0 if save only
param_map('fl_default_aprime') = 0; % if default, next period aprime
param_map('bl_default') = 0; % if borrowing is default allowed

% Assest Grids Toal savings aggregate grid (see discussion on top). 35
% points picked for
% param_map('fl_w_min') = param_map('fl_b_bd'); % but b_bd overrides this
param_map('fl_w_max') = 50;
param_map('it_w_n') = 50;

% Risky Capital Asset Vector see graph below for how it looks graphically
% in principle keep it_k_n the same as it_w_n to have equi-distance points
% in triangle choice grid.
param_map('fl_k_min') = 0;
% param_map('fl_k_max') = (param_map('fl_w_max') - param_map('fl_b_bd'));
param_map('it_ak_n') = param_map('it_w_n'); % grid for a and k the same

% Interpolation fl_coh_interp_grid_gap controls the number of coh points at
% which to solve the model it_c_interp_grid_gap determines the gap between
% consumption terpolation points. For consumption interpolation 10^-4 is
% extremely accurate, there should be no perceptible differences in value
% and policy functions when the it_c_interp_grid_gap <= 0.001 compared to
% actual evaluation
param_map('fl_coh_interp_grid_gap') = 0.1;
% param_map('it_coh_interp_n') = 500;
param_map('it_c_interp_grid_gap') = 10^-4;

% Solution Accuracy
param_map('it_maxiter_val') = 250;
param_map('it_maxiter_dist') = 1000;
param_map('it_trans_power_dist') = 1000;
param_map('st_analytical_stationary_type') = 'eigenvector';
param_map('fl_tol_val') = 10^-5;
param_map('fl_tol_pol') = 10^-5;
param_map('fl_tol_dist') = 10^-5;
param_map('it_tol_pol_nochange') = 25; % number of iterations where policy does not change

%% Setting support_map container

support_map = containers.Map('KeyType','char', 'ValueType','any');

% root directory
[st_root_path] = preamble(false);
st_matimg_path_root = [st_root_path '/m_akz/'];
support_map('st_matimg_path_root') = st_matimg_path_root;

% timer
support_map('bl_time') = true;

% Print Controls
support_map('bl_display_defparam') = false;
support_map('bl_display') = true;
support_map('bl_display_dist') = false;
support_map('it_display_every') = 5; % how often to print results

% Profile Controls
support_map('bl_profile') = false;
support_map('bl_profile_dist') = false; % distribution profile
support_map('st_profile_path') = [st_matimg_path_root '/solve/profile/'];
support_map('st_profile_prefix') = [''];
support_map('st_profile_name_main') = ['_default'];
support_map('st_profile_suffix') = ['_p' num2str(it_subset)];

support_map('bl_post') = false;
% Final Print
support_map('bl_display_final') = false; % print finalized results
support_map('bl_display_final_dist') = false; % print finalized results
support_map('bl_display_final_dist_detail') = false;
support_map('it_display_final_rowmax') = 100; % max row to print (states/iters)
support_map('it_display_final_colmax') = 12; % max col to print (shocks)
it_display_summmat_rowmax = 7;
it_display_summmat_colmax = 7;
support_map('it_display_summmat_rowmax') = it_display_summmat_rowmax;
support_map('it_display_summmat_colmax') = it_display_summmat_colmax;

% Mat File Controls
support_map('bl_mat') = false;
support_map('st_mat_path') = [st_matimg_path_root '/solve/mat/'];
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

% Image Saving Controls (given graphing)
support_map('st_title_prefix') = '';
support_map('bl_img_save') = false;
support_map('st_img_path') = [st_matimg_path_root '/solve/img/'];
support_map('st_img_prefix') = [''];
support_map('st_img_name_main') = ['_default'];
support_map('st_img_suffix') = ['_p' num2str(it_subset) '.png'];

% Sub-function graphing controls
support_map('bl_graph_funcgrids') = false;
support_map('bl_display_funcgrids') = false;
support_map('bl_graph_evf') = false;
support_map('bl_display_evf') = false;

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
        param_map('it_w_n') = 20;
        param_map('it_k_n') = param_map('it_w_n');
        param_map('it_z_n') = 3;
        param_map('it_maxiter_val') = 50;
        param_map('fl_coh_interp_grid_gap') = 0.25;
        param_map('it_c_interp_grid_gap') = 0.001;
        param_map('it_tol_pol_nochange') = 1000;
        support_map('bl_display') = true;
        support_map('it_display_every') = 1;
    end
    if (ismember(it_subset, [2, 4]))
        % close figures
        close all;
        % Main Run
        support_map('bl_time') = true;
        support_map('bl_display_defparam') = true;
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
        param_map('it_w_n') = 20;
        param_map('it_k_n') = param_map('it_w_n');
        param_map('it_z_n') = 3;
        param_map('it_maxiter_val') = 50;
        param_map('fl_coh_interp_grid_gap') = 0.25;
        param_map('it_c_interp_grid_gap') = 0.001;
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
        support_map('bl_display_defparam') = true;
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

        if (ismember(it_subset, [8, 9]))
            support_map('bl_time') = false;
            support_map('bl_display') = false;
            support_map('bl_display_dist') = false;
            support_map('bl_display_final_dist_detail') = true;
            support_map('bl_graph_onebyones') = true;
            support_map('bl_img_save') = false;
            if (ismember(it_subset, [9]))
                % quietly turn off all graphs, only tables
                support_map('bl_display_final_dist_detail') = false;
                support_map('bl_display_defparam') = false;
                support_map('bl_graph_coh_t_coh') = false;
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
    fft_container_map_display(param_map, it_display_summmat_rowmax, it_display_summmat_colmax);
    fft_container_map_display(support_map, it_display_summmat_rowmax, it_display_summmat_colmax);
end

end
