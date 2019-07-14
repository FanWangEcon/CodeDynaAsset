%% Set Model Parameters (Interpolated + Percentage + Risky + Safe Asset + FIBS + RShock)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(varargin)
%% FFS_IPWKBZ_FIBS_SET_DEFAULT_PARAM setting model default parameters

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

%% 1a. Initiate Param_map

param_map = containers.Map('KeyType','char', 'ValueType','any');

% model name
param_map('st_model') = 'ipwkbzr_fibs';


%% 2a. Set Borrowing Control Parameters
% Borrowing Setting 1: Default Allowed, Bridge True, bl_rollover does not matter
% Borrowing Setting 2: Default Allowed, Bridge False, bl_rollover matter
param_map('bl_default') = true; % if borrowing is default allowed
param_map('bl_bridge') = true;
param_map('bl_rollover') = true;
param_map('bl_b_is_principle') = true;

%% 3. Set Interest Rates non-Shock Parameters

% formal informal parameters
% fl_for_br_block are the formal borrowing grid block sizes.
param_map('fl_r_fsv') = 0.025;
param_map('fl_r_fbr') = 0.065;
% see: ffs_for_br_block.m
param_map('st_forbrblk_type') = 'seg3';
param_map('fl_forbrblk_brmost') = -19;
param_map('fl_forbrblk_brleast') = -1;
param_map('fl_forbrblk_gap') = -1.5;

%% 3a. Set Shock 1 Borrowing Interest Rate Parameters
% See
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_gen_discrete_var.html
% fft_gen_discrete_var> for how these parameters will be used to generate a
% discrete random variable for the interest rate.

param_map('st_z_r_infbr_drv_ele_type') = 'unif';
param_map('st_z_r_infbr_drv_prb_type') = 'poiss';
param_map('fl_z_r_infbr_poiss_mean') = 20;
param_map('fl_z_r_infbr_max') = 0.095;
param_map('fl_z_r_infbr_min') = 0.025;
param_map('fl_z_r_infbr_n') = 5;
% param_map('fl_z_r_infbr_max') = 0.095;
% param_map('fl_z_r_infbr_min') = 0.095;
% param_map('fl_z_r_infbr_n') = 1;

%% 3b. Set Shock 2 Productivity Shock Parameters

% Production Function
% Productivity Shock Parameters
param_map('it_z_wage_n') = 15;
param_map('fl_z_wage_mu') = 0;
param_map('fl_z_wage_rho') = 0.8;
param_map('fl_z_wage_sig') = 0.2;

%% 3c. Set Overall Shock Grid Count

param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');

%% Setting support_map container

support_map = containers.Map('KeyType','char', 'ValueType','any');

% Additional Controls
support_map('bl_display_minccost') = false;
support_map('bl_display_infbridge') = false;
support_map('bl_graph_funcgrids') = false;
support_map('bl_graph_funcgrids_detail') = false;
support_map('bl_display_funcgrids') = false;

support_map('bl_graph_forinf_discrete') = true;
support_map('bl_graph_forinf_pol_lvl') = true;
support_map('bl_graph_forinf_pol_pct') = true;

% root directory
[st_root_path] = preamble(false);
st_matimg_path_root = [st_root_path '/m_fibs/'];
support_map('st_matimg_path_root') = st_matimg_path_root;
support_map('st_profile_path') = [st_matimg_path_root '/m_ipwkbzr_solve/profile/'];
support_map('st_mat_path') = [st_matimg_path_root '/m_ipwkbzr_solve/mat/'];
support_map('st_img_path') = [st_matimg_path_root '/m_ipwkbzr_solve/img/'];

%% Display New Parameters

if (bl_display_defparam)

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Display Parameters Specific to IPWKBZR_FIBS')
    disp('it_coh_bridge_perc_n ADDED ON NEXT')
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');

    fft_container_map_display(param_map);
    fft_container_map_display(support_map);
end

%% 3. Merge Parameters Import

[param_map_ipwkbzr_fibs, support_map_ipwkbzr_fibs] = ffs_ipwkbz_set_default_param(it_subset);

% Remove Keys not Relevant for the Interest Rate Shock Model
cl_st_ipwkbzr_keysdrop = {'fl_z_rho', 'fl_z_mu', 'fl_z_sig', ...
                          'fl_r_borr', 'fl_r_save'};
remove(param_map_ipwkbzr_fibs, cl_st_ipwkbzr_keysdrop);

% Merge
param_map = [param_map_ipwkbzr_fibs; param_map];
support_map = [support_map_ipwkbzr_fibs; support_map];

%% Add on based on existing

% Percentage of w that is not for bridge loan, when param_map('bl_bridge') = false
% ar_coh_bridge_perc = [1]
param_map('it_coh_bridge_perc_n') = param_map('it_w_perc_n');

%% Subset Options Adjustments

% close all
close all;

if (ismember(it_subset, [3]))
    % Profile run
elseif (ismember(it_subset, [1,2,4]))
    % Main Run
    if (ismember(it_subset, [1]))

        param_map('it_z_wage_n') = 5;
        param_map('fl_z_r_infbr_n') = 3;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');

        param_map('it_coh_bridge_perc_n') = param_map('it_w_perc_n');

        support_map('bl_graph') = true;
    end
    if (ismember(it_subset, [4]))
    end
end

%% Subset Options for Distribution solutions

if (ismember(it_subset, [7]))
    % Profile run
elseif (ismember(it_subset, [5,6,8,9]))
    % Main Run
        if (it_subset == 5)

        param_map('it_z_wage_n') = 5;
        param_map('fl_z_r_infbr_n') = 3;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');

        param_map('it_coh_bridge_perc_n') = param_map('it_w_perc_n');

        support_map('bl_graph') = true;
    end
    if (ismember(it_subset, [8, 9]))
        if (ismember(it_subset, [9]))
        end

    end

end

%% Display All Parameters

if (bl_display_defparam)

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Display All Parameters with IPWKBZR overriding IPWKBZR')
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');

    fft_container_map_display(param_map);
    fft_container_map_display(support_map);
end

end
