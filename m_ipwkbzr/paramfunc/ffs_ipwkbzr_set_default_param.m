%% Set Model Parameters (Interpolated + Percentage + Risky + Safe Asset + Save + Borrow + R Shock)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [param_map, support_map] = ffs_ipwkbzr_set_default_param(varargin)
%% FFS_IPKBZR_SET_DEFAULT_PARAM setting model default parameters
% Define model parameters, similar to
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_default_param.html
% ffs_akz_set_default_param> see that file for descriptions.
%
% Several changes here: 1, inclusion of percentage based choice grids
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
%   [param_map, support_map] = ffs_ipwkbzr_set_default_param(it_param_set);
%
% @include
%
% *<https://fanwangecon.github.io/CodeDynaAsset/m_ipwkz/paramfunc/html/ffs_ipwkz_set_default_param.html ffs_ipwkz_set_default_param>
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

%% 1. Initiate Param_map

param_map = containers.Map('KeyType','char', 'ValueType','any');

% model name
param_map('st_model') = 'ipwkbzr';

% v(coh, z) interpolation method
param_map('st_v_coh_z_interp_method') = 'method_cell';

%% 2a. Borrowing Default Parameters

param_map('fl_b_bd') = -20;
param_map('fl_w_min') = param_map('fl_b_bd');
param_map('fl_w_max') = 50;
param_map('fl_k_max') = (param_map('fl_w_max') - param_map('fl_b_bd'));

param_map('fl_c_min') = 0.02;
param_map('fl_default_wprime') = 0; % wprime not a prime
param_map('bl_default') = true; % if borrowing is default allowed

%% 2b. Set Shock 1 Borrowing Interest Rate Parameters
% See
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_gen_discrete_var.html
% fft_gen_discrete_var> for how these parameters will be used to generate a
% discrete random variable for the interest rate.

param_map('st_z_r_borr_drv_ele_type') = 'unif';
param_map('st_z_r_borr_drv_prb_type') = 'poiss';
param_map('fl_z_r_borr_poiss_mean') = 20;
param_map('fl_z_r_borr_max') = 0.095;
param_map('fl_z_r_borr_min') = 0.025;
param_map('fl_z_r_borr_n') = 5;
% param_map('fl_z_r_borr_max') = 0.095;
% param_map('fl_z_r_borr_min') = 0.095;
% param_map('fl_z_r_borr_n') = 1;

%% 2c. Set Shock 2 Productivity Shock Parameters

% Production Function
% Productivity Shock Parameters
param_map('it_z_wage_n') = 11;
param_map('fl_z_wage_mu') = 0;
param_map('fl_z_wage_rho') = 0.8;
param_map('fl_z_wage_sig') = 0.2;

%% 2d. Set Overall Shock Grid Count

param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

%% Setting support_map container

support_map = containers.Map('KeyType','char', 'ValueType','any');

% root directory
[st_root_path] = preamble(false);
st_matimg_path_root = [st_root_path '/m_ipwkbzr/'];
support_map('st_matimg_path_root') = st_matimg_path_root;
support_map('st_profile_path') = [st_matimg_path_root '/solve/profile/'];
support_map('st_mat_path') = [st_matimg_path_root '/solve/mat/'];
support_map('st_img_path') = [st_matimg_path_root '/solve/img/'];
support_map('st_mat_test_path') = [st_matimg_path_root '/test/ff_ipwkbzr_ds_vecsv/mat/'];

%% Subset Options
%
% # it_subset = 1 is basic invoke quick test
% # it_subset = 2 is main invoke
% # it_subset = 3 is profiling invoke
% # it_subset = 4 is matlab publish.
%

% close figures
close all;

if (ismember(it_subset, [3]))
    % Profile run
elseif (ismember(it_subset, [1,2,4]))
    % Main Run
    if (ismember(it_subset, [1]))
        % TEST quick
        param_map('it_z_wage_n') = 5;
        param_map('fl_z_r_borr_n') = 3;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
    end
    if (ismember(it_subset, [4]))
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

if (ismember(it_subset, [7]))
    % Profile run
elseif (ismember(it_subset, [5,6,8,9]))
    % Main Run
    if (ismember(it_subset, [5]))
        % TEST quick (need to enough to have distribution)
        param_map('it_z_wage_n') = 5;
        param_map('fl_z_r_borr_n') = 3;
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
    end

    if (ismember(it_subset, [8, 9]))
        if (ismember(it_subset, [9]))
            % quietly turn off all graphs, only tables
        end
    end
end

%% Display New Parameters

if (bl_display_defparam)

    disp('----------------------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('Display Parameters Specific to IPWKBZR')
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');

    fft_container_map_display(param_map);
    fft_container_map_display(support_map);
end

%% 3. Merge Parameters Import

[param_map_ipwkbz, support_map_ipwkbz] = ffs_ipwkz_set_default_param(it_subset);

% Remove Keys not Relevant for the Interest Rate Shock Model
cl_st_ipwkbz_keysdrop = {'fl_z_rho', 'fl_z_mu', 'fl_z_sig', ...
                         'fl_r_borr'};
remove(param_map_ipwkbz, cl_st_ipwkbz_keysdrop);

% Merge
param_map = [param_map_ipwkbz; param_map];
support_map = [support_map_ipwkbz; support_map];

% if (bl_display_defparam)
%
%     disp('----------------------------------------');
%     disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
%     disp('Display Parameters imported from IPWKBZ')
%     disp('Dropped Keys:')
%     disp(cl_st_ipwkbz_keysdrop)
%     disp('xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx');
%
%     fft_container_map_display(param_map_ipwkbz);
%     fft_container_map_display(support_map_ipwkbz);
% end


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
