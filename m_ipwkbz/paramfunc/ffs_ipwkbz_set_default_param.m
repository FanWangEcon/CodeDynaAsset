%% Set Model Parameters (Interpolated + Percentage + Risky + Safe Asset + Save + Borrow)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [param_map, support_map] = ffs_ipwkbz_set_default_param(varargin)
%% FFS_IPKBZ_SET_DEFAULT_PARAM setting model default parameters
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
%   [param_map, support_map] = ffs_ipwkbz_set_default_param(it_param_set);
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
param_map('st_model') = 'ipwkbz';

% v(coh, z) interpolation method
param_map('st_v_coh_z_interp_method') = 'method_cell';

%% 2. Borrowing Default Parameters
param_map('fl_b_bd') = -20;
param_map('fl_c_min') = 0.02;
param_map('fl_default_wprime') = 0; % wprime not a prime
param_map('bl_default') = true; % if borrowing is default allowed

%% Setting support_map container

support_map = containers.Map('KeyType','char', 'ValueType','any');

% root directory
[st_root_path] = preamble(false);
st_matimg_path_root = [st_root_path '/m_ipwkbzr/'];
support_map('st_matimg_path_root') = st_matimg_path_root;
support_map('st_profile_path') = [st_matimg_path_root '/solve/profile/'];
support_map('st_mat_path') = [st_matimg_path_root '/solve/mat/'];
support_map('st_img_path') = [st_matimg_path_root '/solve/img/'];
support_map('st_mat_test_path') = [st_matimg_path_root '/test/ff_ipwkbz_ds_vecsv/mat/'];

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
