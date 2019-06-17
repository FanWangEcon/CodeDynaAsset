%% Test Borrowing Grid Default (Save + Borr Dynamic Programming Problem)
% Testing the
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
% program which generates the borrowing and savings choice grid required
% for solving the borrowing and savings dynamic programming problem:
% <https://fanwangecon.github.io/CodeDynaAsset/m_abz/solve/html/ff_abz_vf_vecsv.html
% ff_abz_vf_vecsv>. 
%
% Below I show the choice grid for (1) savings only (2) borrowing without
% default, which means cash-on-hand next period must be always positive (3)
% borrowing grid with default, which extends up to negative asset levels
% where there is at least one state of high enough shock in which the
% borrower can still repay debt without defaulting. 
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

close all

% Production Function
fl_a_max = 2;
it_a_n = 50;
fl_r_save = 0.10; % 10 percent savings interest rate
fl_r_borr = 0.50;
fl_w = 3;

% Display
bl_graph_funcgrids = true;
bl_display_funcgrids = false;

%% Generate Savings a grid

% Not default parameters, but parameters that generate defaults
it_param_set = 4;
bl_input_override = true;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% graph and display
support_map('bl_graph_funcgrids') = bl_graph_funcgrids;
support_map('bl_display_funcgrids') = bl_display_funcgrids;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% control saving, borrowing, default
param_map('fl_b_bd') = 0;
param_map('bl_default') = 0;

% run program
[armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);


%% Generate Borrowing A Grid without Default

% Not default parameters, but parameters that generate defaults
it_param_set = 4;
bl_input_override = true;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% graph and display
support_map('bl_graph_funcgrids') = bl_graph_funcgrids;
support_map('bl_display_funcgrids') = bl_display_funcgrids;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% control saving, borrowing, default
param_map('fl_b_bd') = -20;
param_map('bl_default') = 0;

% run program
[armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);


%% Generate Borrowing A Grid with Default

% Not default parameters, but parameters that generate defaults
it_param_set = 4;
bl_input_override = true;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% graph and display
support_map('bl_graph_funcgrids') = bl_graph_funcgrids;
support_map('bl_display_funcgrids') = bl_display_funcgrids;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% control saving, borrowing, default
param_map('fl_b_bd') = -20;
param_map('bl_default') = 1;

% run program
[armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);


%% Generate Borrowing A Grid with Binding Exo Borrowing

% Not default parameters, but parameters that generate defaults
it_param_set = 4;
bl_input_override = true;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);

% graph and display
support_map('bl_graph_funcgrids') = bl_graph_funcgrids;
support_map('bl_display_funcgrids') = bl_display_funcgrids;

% shared parameters
param_map('fl_a_min') = 0;
param_map('fl_a_max') = fl_a_max;
param_map('it_a_n') = it_a_n;
param_map('fl_r_save') = fl_r_save;
param_map('fl_r_borr') = fl_r_borr;
param_map('fl_w') = fl_w;

% control saving, borrowing, default
param_map('fl_b_bd') = -8;
param_map('bl_default') = 1;

% run program
[armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map, bl_input_override);
