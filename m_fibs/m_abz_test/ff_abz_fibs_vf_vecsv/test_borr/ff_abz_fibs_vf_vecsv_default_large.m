%% Test Borrowing with Default (Standard Grid, For + Inf + Borr + Save DP)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Testing the
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_fibs_vf_vecsv.html ff_abz_fibs_vf_vecsv>
% program for solving the savings and borrowing dynamic programming problem.
%
% Test the model by adjusting (1) borrowing bound (2) if default is allowed
% (3) when default is allowed adjusting the level of minimum consumption
% given default.
%
% Here graphical results are printed out for large grid solution.
%
% @seealso
%
% * _BORROW_: borrow and default small grid testing: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_test/ff_abz_fibs_vf_vecsv/test_borr/html/ff_abz_fibs_vf_vecsv_default_small.html ff_abz_fibs_vf_vecsv_default_small>
% * _BORROW_: borrow and default large grid testing: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_test/ff_abz_fibs_vf_vecsv/test_borr/html/ff_abz_fibs_vf_vecsv_default_large.html ff_abz_fibs_vf_vecsv_default_large>
% * _PRECISION_: borr + save quick vs benchmark testing: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_test/ff_abz_fibs_vf_vecsv/test_precision/html/fsi_abz_fibs_vf_vecsv_main.html fsi_abz_fibs_vf_vecsv_main>
% * _PRECISION_: borr + save asset grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_test/ff_abz_fibs_vf_vecsv/test_precision/html/fsi_abz_fibs_vf_vecsv_a_n.html fsi_abz_fibs_vf_vecsv_a_n>
% * _PRECISION_: borr + save shock grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_test/ff_abz_fibs_vf_vecsv/test_precision/html/fsi_abz_fibs_vf_vecsv_z_n.html fsi_abz_fibs_vf_vecsv_z_n>
%

%% Set Shared Parameters

close all;

it_param_set = 4;

% Shared parameters
fl_a_max = 50;
it_a_n = 750;
fl_r_save = 0.025; % 10 percent savings interest rate
fl_r_borr = 0.035;
fl_w = 1;
it_maxiter_val = 1000;
it_z_n = 15;

% Display Etc.
bl_display = false;
bl_post = true;
bl_display_final = true;

%% Simulate Savinags only
% default not allowed and borrowing bound = 0

[param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);

% control saving, borrowing, default
param_map('fl_b_bd') = 0;
param_map('bl_default') = false;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% solution
param_map('it_z_n') = it_z_n;
param_map('it_maxiter_val') = it_maxiter_val;

% Display Parameters
support_map('bl_display') = bl_display;
support_map('bl_post') = bl_post;
support_map('bl_display_final') = bl_display_final;

% Call Program
ff_abz_fibs_vf_vecsv(param_map, support_map);

%% Simulate Save/Borrow No Default
% Default not allowed, but borrowing is allowed
%

[param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);

% allow for borrowing, when using a very large negative value for fl_b_bd
% that means the natural borrowing constraint will bind.
param_map('fl_b_bd') = -100;
param_map('bl_default') = false;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% solution
param_map('it_z_n') = it_z_n;
param_map('it_maxiter_val') = it_maxiter_val;

% Display Parameters
support_map('bl_display') = bl_display;
support_map('bl_post') = bl_post;
support_map('bl_display_final') = bl_display_final;

% Call Program
ff_abz_fibs_vf_vecsv(param_map, support_map);

% Snap
snapnow;

close all;

%% Save/Borrow, Can Default, No Bridge: cmin = 0.00001, choose not to default
% Default is allowed, borrowing is allowed, very low minimum consumption
% value when default happens.
%

[param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);

% allow for borrowing, when using a very large negative value for fl_b_bd
% that means the natural borrowing constraint will bind.
param_map('fl_b_bd') = -20;
param_map('bl_default') = true;
param_map('fl_c_min') = 0.00001;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% solution
param_map('it_z_n') = it_z_n;
param_map('it_maxiter_val') = it_maxiter_val;

% Display Parameters
support_map('bl_display') = bl_display;
support_map('bl_post') = bl_post;
support_map('bl_display_final') = bl_display_final;

% Call Program
ff_abz_fibs_vf_vecsv(param_map, support_map);

% Snap
snapnow;

close all;

%% Save/Borrow, Can Default, No Bridge: cmin = 1, degenerate, borrow to max
% Default is allowed, borrowing is allowed, minimum consumption = 1, very
% high. Always defaulting, distribution degenerate.
%

[param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);

% allow for borrowing, when using a very large negative value for fl_b_bd
% that means the natural borrowing constraint will bind.
param_map('fl_b_bd') = -20;
param_map('bl_default') = 1;
param_map('fl_c_min') = 1;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% solution
param_map('it_z_n') = it_z_n;
param_map('it_maxiter_val') = it_maxiter_val;

% Display Parameters
support_map('bl_display') = bl_display;
support_map('bl_post') = bl_post;
support_map('bl_display_final') = bl_display_final;


% Call Program
ff_abz_fibs_vf_vecsv(param_map, support_map);

% Snap
snapnow;

close all;

%% Save/Borrow, Can Default, No Bridge: cmin = 0.001, default takes place
% Default is allowed, borrowing is allowed, minimum consumption = 0.01,
% higher than previous.
%

[param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);

% allow for borrowing, when using a very large negative value for fl_b_bd
% that means the natural borrowing constraint will bind.
param_map('fl_b_bd') = -20;
param_map('bl_default') = 1;
param_map('fl_c_min') = 0.001;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% solution
param_map('it_z_n') = it_z_n;
param_map('it_maxiter_val') = it_maxiter_val;

% Display Parameters
support_map('bl_display') = bl_display;
support_map('bl_post') = bl_post;
support_map('bl_display_final') = bl_display_final;


% Call Program
ff_abz_fibs_vf_vecsv(param_map, support_map);

% Snap
snapnow;

close all;

%% Save/Borrow, Can Default, Yes Bridge: cmin = 0.001, default takes place
% Default is allowed, borrowing is allowed, minimum consumption = 0.01,
% higher than previous.
%

[param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);

% allow for borrowing, when using a very large negative value for fl_b_bd
% that means the natural borrowing constraint will bind.
param_map('fl_b_bd') = -20;
param_map('bl_default') = 1;
param_map('bl_bridge') = true;
param_map('fl_c_min') = 0.001;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% solution
param_map('it_z_n') = it_z_n;
param_map('it_maxiter_val') = it_maxiter_val;

% Display Parameters
support_map('bl_display') = bl_display;
support_map('bl_post') = bl_post;
support_map('bl_display_final') = bl_display_final;


% Call Program
ff_abz_fibs_vf_vecsv(param_map, support_map);

% Snap
snapnow;

close all;

%% Save/Borrow, Can Default, No Bridge, No Rollover: cmin = 0.001, default takes place
% Default is allowed, borrowing is allowed, minimum consumption = 0.01,
% higher than previous.
%

[param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);

% allow for borrowing, when using a very large negative value for fl_b_bd
% that means the natural borrowing constraint will bind.
param_map('fl_b_bd') = -20;
param_map('bl_default') = 1;
param_map('bl_bridge') = false;
param_map('bl_rollover') = false;
param_map('fl_c_min') = 0.001;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% solution
param_map('it_z_n') = it_z_n;
param_map('it_maxiter_val') = it_maxiter_val;

% Display Parameters
support_map('bl_display') = bl_display;
support_map('bl_post') = bl_post;
support_map('bl_display_final') = bl_display_final;


% Call Program
ff_abz_fibs_vf_vecsv(param_map, support_map);

% Snap
snapnow;

close all;
clear all;

