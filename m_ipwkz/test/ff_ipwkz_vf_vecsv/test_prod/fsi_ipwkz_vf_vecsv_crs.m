%% IPWKZ Stock vs Risky Capital Investment

close all

%% Simulate Risky Asset with CRS and full depreciation, IID shocks
% CRRA utility with Constant return to scale risky asset and safe asset at
% the parameters below, constant share of asset in riksy asset, at the
% parameters below, shock does not have additional effects except through
% coh so plotting conditional on cash-on-hand, there is no difference
% across shocks, hence solve at only it_z_n == 3. Note we have to have many
% grid points for it_w_n to get the exact solution. For pratical purposes,
% this does not need to be as high.

it_param_set = 4;
[param_map, support_map] = ffs_ipwkz_set_default_param(it_param_set);

% Simulation Accuracy
param_map('it_w_perc_n') = 500;
param_map('it_ak_perc_n') = param_map('it_w_perc_n');
param_map('it_z_n') = 15;

param_map('fl_coh_interp_grid_gap') = 0.01;
param_map('fl_w_interp_grid_gap') = 0.01;
param_map('it_c_interp_grid_gap') = 10^-4;

% Iteration Limits
% param_map('it_maxiter_val') = 2;

% Turn on 2nd stage graphs
support_map('bl_graph_evf') = false;
support_map('bl_display_evf') = false;

% Production Function Parameters
% note shock is log normal
param_map('fl_Amean') = 1.0265;
param_map('fl_alpha') = 1;
param_map('fl_delta') = 1;
param_map('fl_r') = 0.03;
param_map('fl_w') = 0.20;

% Shock Parameter, iid shocks
param_map('fl_z_rho') = 0;
param_map('fl_z_sig') = 0.05;

% Display Parameters
support_map('bl_display') = false;
support_map('bl_display_final') = true;
support_map('bl_time') = true;
% support_map('bl_profile') = false;


% Call Program
ff_ipwkz_vf_vecsv(param_map, support_map);

%% Simulate Risky Asset with CRS and full depreciation, Persistent Shocks
% Now Shocks matter conditional on coh(z) still.

it_param_set = 4;
[param_map, support_map] = ffs_ipwkz_set_default_param(it_param_set);

% Simulation Accuracy
param_map('it_w_perc_n') = 500;
param_map('it_ak_perc_n') = param_map('it_w_perc_n');
param_map('it_z_n') = 15;

param_map('fl_coh_interp_grid_gap') = 0.01;
param_map('fl_w_interp_grid_gap') = 0.01;
param_map('it_c_interp_grid_gap') = 10^-4;

% Production Function Parameters
% note shock is log normal
param_map('fl_Amean') = 1.0265;
param_map('fl_alpha') = 1;
param_map('fl_delta') = 1;
param_map('fl_r') = 0.03;
param_map('fl_w') = 0.20;

% Shock Parameter, iid shocks
param_map('fl_z_rho') = 0.05;
param_map('fl_z_sig') = 0.05;

% Display Parameters
support_map('bl_display') = false;
support_map('bl_display_final') = false;
support_map('bl_time') = true;
% support_map('bl_profile') = false;

% Call Program
ff_ipwkz_vf_vecsv(param_map, support_map);
