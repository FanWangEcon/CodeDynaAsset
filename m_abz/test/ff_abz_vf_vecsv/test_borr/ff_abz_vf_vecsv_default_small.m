%% Test Borrowing with Default (Small Grid, Save + Borr Dynamic Programming Problem)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Testing the
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf_vecsv.html ff_abz_vf_vecsv>
% program for solving the savings and borrowing dynamic programming problem.
%
% Test the model by adjusting (1) borrowing bound (2) if default is allowed
% (3) when default is allowed adjusting the level of minimum consumption
% given default.
%
% Here tabular solutions are printed out for small grid solution.
%
% @seealso
%
% * _BORROW GRID_: borrowing choice grid with default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ffs_abz_get_funcgrid/test_borr/html/ffs_abz_get_funcgrid_defnodfalt.html ffs_abz_get_funcgrid_defnodfalt>
% * _BORROW_: borrow and default small grid testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_borr/html/ff_abz_vf_vecsv_default_small.html ff_abz_vf_vecsv_default_small>
% * _BORROW_: borrow and default large grid testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_borr/html/ff_abz_vf_vecsv_default_large.html ff_abz_vf_vecsv_default_large>
% * _PRECISION_: borr + save quick vs benchmark testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_precision/html/fsi_abz_vf_vecsv_main.html fsi_abz_vf_vecsv_main>
% * _PRECISION_: borr + save asset grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_precision/html/fsi_abz_vf_vecsv_a_n.html fsi_abz_vf_vecsv_a_n>
% * _PRECISION_: borr + save shock grid count testing: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_abz_vf_vecsv/test_precision/html/fsi_abz_vf_vecsv_z_n.html fsi_abz_vf_vecsv_z_n>
%

%% Set Shared Parameters

close all;

it_param_set = 1;

% Shared parameters
fl_a_max = 50;
it_a_n = 60;
fl_r_save = 0.025; % 10 percent savings interest rate
fl_r_borr = 0.035;
fl_w = 1;
it_maxiter_val = 1000;
it_z_n = 7;

% Display Etc.
bl_display = false;
bl_post = true;
bl_display_final = true;

%% Simulate Savinags only
% for the output optimal choice matrix, look at each row, if the next
% period a at the lowest level of shock is below current a, that means
% there is nonzero chance of going to lower asset level from this point.
% Similarly for higher point. If at a current asset point, for all shock
% levels, next period a choice is equal or below current period asset
% level, this asset level will have zero probability in the stationary
% distribution. For the distribution to be non-degenerate, there must be no
% rows where a'(a,z) == a for all z. see
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ffs_abz_get_funcgrid/test_borr/html/ffs_abz_get_funcgrid_nodefault.html
% ffs_abz_get_funcgrid_nodefault> for additional information.

[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% control saving, borrowing, default
param_map('fl_b_bd') = 0;
param_map('bl_default') = 0;

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
ff_abz_vf_vecsv(param_map, support_map);

%% Simulate Save/Borrow No Default
% See each asset state row, if for negative asset levels, at all levels of
% shocks, a'(a,z) >= a for all a <0 for all z, that means there will be no
% borrowing in the stationary distribution. Similarly if a'(a,z) <= a for
% all a > 0 and all z, there will be no savings in the stationary
% distribution. For the distribution to be non-degenerate, there must be no
% rows where a'(a,z) == a for all z.
%

[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% allow for borrowing, when using a very large negative value for fl_b_bd
% that means the natural borrowing constraint will bind.
param_map('fl_b_bd') = -100;
param_map('bl_default') = 0;

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
ff_abz_vf_vecsv(param_map, support_map);

%% Simulate Save/Borrow, Can Default: cmin = 0.00001, choose not to default
% Allow for default, but very low cmin.
%
% see that the resulting policy function does not allow distribution to
% exceed below -9.32203, the additional borrowing allowed for default does
% not matter because utility from default is so terrible with c min =
% 0.00001.
%

[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% allow for borrowing, when using a very large negative value for fl_b_bd
% that means the natural borrowing constraint will bind.
param_map('fl_b_bd') = -20;
param_map('bl_default') = 1;
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
ff_abz_vf_vecsv(param_map, support_map);

%% Simulate Save/Borrow, Can Default: cmin = 0.001, default takes place
% now borrowing level feasible increases significantly given minimum
% consumption and default. see
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ffs_abz_get_funcgrid/test_borr/html/ffs_abz_get_funcgrid_nodefault.html
% ffs_abz_get_funcgrid_nodefault> for additional information.
%
% higher cmin than before, see that the lowest asset level reached is now
% lower. See that around -10, there is no a' below the 45 degree line in
% the asset choice graph.
%

[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

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
ff_abz_vf_vecsv(param_map, support_map);


%% Simulate Save/Borrow, Can Default: cmin = 1, degenerate, borrow to max
% with cmin = 1, we have garanteed consumption regardless debt owed. Now
% the distribution generated by the policy functions will be degenerate. In
% the negative asset region, all choices move towards the defaulting
% states.
%

[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

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
ff_abz_vf_vecsv(param_map, support_map);
