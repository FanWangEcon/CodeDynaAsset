%% Test FFS_ABZ_GET_FUNCGRID handling of defaults
% Here I test
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_abz_get_funcgrid.html
% ffs_abz_get_funcgrid> the choice grid changes when we have savings vs
% borrowing without defaults vs borrowing with default.

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
