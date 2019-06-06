%% IPWKZ_EVF CRS investment Percentage based grid vs Grid
% One type of risky investment is the stocks kind, with constant returns to
% scale and full depreciation. We have the
% *ipwkz*=interpolated-percentage-grid-based-2-stage-solution. What is the
% point of this algorithm, vs, *iwkz* which seems pretty good?. When we
% solve the risky physical investment choice, the 2nd stage of the two
% problems are not very different it seems, but when the level of optimal
% risky investment changes smoothly as aggregate savings increase (which is
% the case when we have risky stocks) it makes a very big difference
% whether we use *ipwkz* or *iwkz*.
%
% Below I call the two 2nd stage solution files when we have stocks and you
% can see clearly the percentage solution is very accurate at low levels of
% aggregate savings choices. The levels choices don't look visually
% different, but when we do percentage, which allows for zooming in to
% optimal choices with low w levels, you see the difference.
%
% Compare to
% <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkz/test/ff_ipwkz_evf/test_prod/html/fsi_ipwkz_vf_vecsv_physicalk.html
% fsi_ipwkz_vf_vecsv_physicalk> which is the companion file.
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_wkz_evf.html ff_wkz_evf>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_ipwkz/solve/html/ff_ipwkz_evf.html ff_ipwkz_evf>
%


close all

% Production Function
fl_Amean = 1.0265;
fl_alpha = 1;
fl_delta = 1;
fl_r = 0.03;
fl_w = 0;

% Shock Parameter, iid shocks
fl_z_rho = 0.005;
fl_z_sig = 0.05;

% Choice Min and Max
fl_b_bd = 0;
fl_w_max = 50;

% Grids, level grid is upper triangle, percentage grid is full N by N
it_ak_n = 250;
it_ak_perc_n = round(sqrt(it_ak_n*(it_ak_n-1)/2+it_ak_n));

% Display
bl_graph_evf = true;
bl_display_evf = false;

%% Solve 2nd Stage Percentage Grid k(w,z) choices, ipkwz
% If it is important to get low level choices percentage levels properly,
% one need to use the percentage grid solution. Below. I show the
% percentage grid based solution when risky investment has CRS and full
% depreciation, the "stock" example. This investment. We only solve the
% second stage problem.
%

close all;

% Not default parameters, but parameters that generate defaults
it_param_set = 4;
bl_input_override = true;
[param_map, support_map] = ffs_ipwkz_set_default_param(it_param_set);

support_map('bl_graph_evf') = bl_graph_evf;
support_map('bl_display_evf') = bl_display_evf;

% 177 because 177^2 = 31375 approximately, which is the grid for the level
% grid below, for fair comparison
param_map('fl_b_bd') = fl_b_bd;
param_map('fl_w_max') = fl_w_max;
param_map('fl_w_min') = param_map('fl_b_bd');
param_map('it_ak_perc_n') = it_ak_perc_n;
param_map('fl_w_interp_grid_gap') = (param_map('fl_w_max')-param_map('fl_b_bd'))/param_map('it_ak_perc_n');

% note shock is log normal
param_map('fl_Amean') = fl_Amean;
param_map('fl_alpha') = fl_alpha;
param_map('fl_delta') = fl_delta;
param_map('fl_r') = fl_r;
param_map('fl_w') = fl_w;

% Shock Parameter, iid shocks
param_map('fl_z_rho') = fl_z_rho;
param_map('fl_z_sig') = fl_z_sig;

[armt_map, func_map] = ffs_ipwkz_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override

% Generating Defaults
params_group = values(armt_map, {'ar_a_meshk', 'ar_k_mesha', 'ar_z'});
[ar_a_meshk, ar_k_mesha, ar_z] = params_group{:};
params_group = values(func_map, {'f_util_standin'});
[f_util_standin] = params_group{:};
mt_val = f_util_standin(ar_z, ar_a_meshk, ar_k_mesha);

% Call Program
ff_ipwkz_evf(mt_val, param_map, support_map, armt_map, bl_input_override);

%% Solve 2nd Stage Fixed Level Grid k(w,z) choices, akz/wkz/iwkz
% This is the grid in the akz/wkz/iwkz problems, we have a fixed grid, not
% a percentage based grid.

% Not default parameters, but parameters that generate defaults
it_param_set = 4;
bl_input_override = true;
[param_map, support_map] = ffs_akz_set_default_param(it_param_set);
support_map('bl_graph_evf') = bl_graph_evf;
support_map('bl_display_evf') = bl_display_evf;

param_map('fl_b_bd') = fl_b_bd;
param_map('fl_w_max') = fl_w_max;
param_map('fl_w_min') = param_map('fl_b_bd');
param_map('it_ak_n') = it_ak_n;
param_map('it_w_n') = param_map('it_ak_n');
% this requires 250*(250-1)/2+250 = 31375 solution points

% note shock is log normal
param_map('fl_Amean') = fl_Amean;
param_map('fl_alpha') = fl_alpha;
param_map('fl_delta') = fl_delta;
param_map('fl_r') = fl_r;
param_map('fl_w') = fl_w;

% Shock Parameter, iid shocks
param_map('fl_z_rho') = fl_z_rho;
param_map('fl_z_sig') = fl_z_sig;

[armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override

% Generating Defaults
params_group = values(armt_map, {'ar_a_meshk', 'ar_k_mesha', 'ar_z'});
[ar_a_meshk, ar_k_mesha, ar_z] = params_group{:};
params_group = values(func_map, {'f_util_standin'});
[f_util_standin] = params_group{:};
mt_val = f_util_standin(ar_z, ar_a_meshk, ar_k_mesha);

% Call Program
ff_wkz_evf(mt_val, param_map, support_map, armt_map, bl_input_override);
