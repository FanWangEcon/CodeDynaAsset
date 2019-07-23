%% Set Parameters One Asset Borrowing and Savings
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [param_map, support_map] = ffs_abz_set_default_param(varargin)
%% FFS_ABZ_SET_DEFAULT_PARAM setting model default parameters
% Two groups of default parameters stored in container maps. This function
% initiliazes parameters.
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
%   [param_map, support_map] = ffs_abz_set_default_param(it_param_set);
%
% @seealso
%
% * initialize paramters *az*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html ffs_az_set_default_param>
% * initialize paramters *abz*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_set_default_param.html ffs_abz_set_default_param>
% * initialize functions: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_set_functions.html ffs_abz_set_functions>
% * set asset grid: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_gen_borrsave_grid.html ffs_abz_gen_borrsave_grid>
% * set shock borrow rate: <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_gen_discrete_var.html fft_gen_discrete_var>
% * set shock wage: <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/tools/ffto_gen_tauchen_jhl.m ffto_gen_tauchen_jhl>
% * gateway function processing grid, paramters, functions: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
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

%% 1. Initiate Param_map

param_map = containers.Map('KeyType','char', 'ValueType','any');

% model name
param_map('st_model') = 'abz';

% Preferences
param_map('fl_crra') = 1.5;
param_map('fl_beta') = 0.94;
param_map('fl_nan_replace') = -99999;

%% 2a. Set Borrowing Control Parameters

% Borrowing
% fl_default_aprime is the next period asset level
% households face if they default.
param_map('fl_b_bd') = -20; % borrow bound, = 0 if save only
param_map('fl_default_aprime') = 0; % if default, next period aprime
param_map('bl_default') = true; % if borrowing is default allowed

% is save/borr choice principle or principle + interest, matters for
% borrowing grid generation program. the *abz* problem is written with
% asset choice as principle only, the _abz_fibs_ problems are written as
% priniple + interest as the state, so there this should be false.
param_map('bl_b_is_principle') = true;

% Minimum Consumption, c_min is for default, when c < 0, replace utility
% with fl_nan_replace.
param_map('fl_c_min') = 0.01;

%% 2b. Set Asset Grid Parameters
% see
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_gen_borrsave_grid.html
% ffs_abz_gen_borrsave_grid> for how these borrowing/saving grid parameters
% will be used.

% Savings
param_map('fl_a_min') = 0; % if there is minimum savings requirement
param_map('fl_a_max') = 50;
param_map('bl_loglin') = false; % log lin threshold structure
param_map('fl_loglin_threshold') = 1; % dense points before 1
param_map('it_a_n') = 750;
% Prices
param_map('fl_w') = 1.28;

%% 3a. Set Interest Rates non-Shock Parameters

param_map('fl_r_save') = 0.025;

%% 4a. Set Shock 1 Borrowing Interest Rate Parameters
% See
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_gen_discrete_var.html
% fft_gen_discrete_var> for how these parameters will be used to generate a
% discrete random variable for the interest rate.

param_map('st_z_r_borr_drv_ele_type') = 'unif';
param_map('st_z_r_borr_drv_prb_type') = 'poiss';
param_map('fl_z_r_borr_poiss_mean') = 10;
param_map('fl_z_r_borr_max') = 0.095;
param_map('fl_z_r_borr_min') = 0.025;
param_map('fl_z_r_borr_n') = 5;
% param_map('fl_z_r_borr_max') = 0.095;
% param_map('fl_z_r_borr_min') = 0.095;
% param_map('fl_z_r_borr_n') = 1;

%% 4b. Set Shock 2 Wage Shock Parameters
% See
% <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/tools/ffto_gen_tauchen_jhl.m
% ffto_gen_tauchen_jhl> for standard implementation of ar1 shock process
% using these parameters.

param_map('it_z_wage_n') = 15;
param_map('fl_z_wage_mu') = 0;
param_map('fl_z_wage_rho') = 0.8;
param_map('fl_z_wage_sig') = 0.2;

%% 4c. Set Overall Shock Grid Count

param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

%% 5. Set Solution Control Parameters

% Solution Accuracy
param_map('it_maxiter_val') = 1000;
param_map('it_maxiter_dist') = 1000;
param_map('st_analytical_stationary_type') = 'eigenvector'; % could be eigenvector, projection, power
param_map('it_trans_power_dist') = 1000;
param_map('fl_tol_val') = 10^-5;
param_map('fl_tol_pol') = 10^-5;
param_map('fl_tol_dist') = 10^-5;
param_map('it_tol_pol_nochange') = 25; % number of iterations where policy does not change

%% 6. Setting support_map container

support_map = containers.Map('KeyType','char', 'ValueType','any');
% root directory
[st_root_path] = preamble(false);
st_matimg_path_root = [st_root_path '/m_abz/'];
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
support_map('bl_display_final_dist_detail') = false; % print finalized results
support_map('it_display_final_rowmax') = 100; % max row to print (states/iters)
support_map('it_display_final_colmax') = 15; % max col to print (shocks)
it_display_summmat_rowmax = 5;
it_display_summmat_colmax = 5;
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
% Testing Path
support_map('st_mat_test_path') = [st_matimg_path_root '/test/ff_az_ds_vecsv/mat/'];

% Sub-function graphing controls
support_map('bl_graph_funcgrids') = false;
support_map('bl_display_funcgrids') = false;
support_map('bl_display_defparam') = false;

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
        param_map('it_z_wage_n') = 3;
        param_map('fl_z_r_borr_n') = 2;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
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
            support_map('bl_display') = false;
            support_map('bl_graph_onebyones') = true;
            support_map('bl_img_save') = false;
        end
    end
    if (ismember(it_subset, [3]))
        % Profile run
        support_map('bl_profile') = true;
        support_map('bl_display') = false; % don't print
        support_map('bl_display_dist') = false; % don't print
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

        param_map('it_z_wage_n') = 5;
        param_map('fl_z_r_borr_n') = 2;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

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

%% Display
if (bl_display_defparam)
    fft_container_map_display(param_map, it_display_summmat_rowmax, it_display_summmat_colmax);
    fft_container_map_display(support_map, it_display_summmat_rowmax, it_display_summmat_colmax);
end

end
