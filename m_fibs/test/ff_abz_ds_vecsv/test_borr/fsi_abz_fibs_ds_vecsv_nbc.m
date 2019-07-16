%% For + Inf + Borr + Save, No Default Testings
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% Testing borrowing and savings cases when default is not allowed
%
% In principle, changes in cmin should not impact results since default
% does not occur here. cmin should not enter into the solution at all. 
%
% @seealso
%

%% Common Parameters
bl_default = false;
fl_r_inf_low = 0.07;
fl_r_inf_high = 0.10;
fl_r_fbr = 0.06;

%% Set Parameters where formal and informal interest rates are the same
param_map_group_a = containers.Map('KeyType','char', 'ValueType','any');

% Key Borrowing Controls
param_map_group_a('bl_default') = bl_default;

% Interset Rates
param_map_group_a('fl_r_inf') = fl_r_inf_low;
param_map_group_a('fl_r_inf_bridge') = fl_r_inf_low;
param_map_group_a('fl_r_fbr') = fl_r_fbr;

%% Set Parameters where formal and informal interest rates are the different
param_map_group_b = containers.Map('KeyType','char', 'ValueType','any');

% Key Borrowing Controls
param_map_group_b('bl_default') = bl_default;

% Interset Rates
param_map_group_b('fl_r_inf') = fl_r_inf_high;
param_map_group_b('fl_r_inf_bridge') = fl_r_inf_high;
param_map_group_b('fl_r_fbr') = fl_r_fbr;

%% *Case A1* No Default, Formal Menu (formal can rollover), fl_r_inf (0.07) > fl_r_fbr (0.6)
% Same Formal and informal rate

param_map_group_a('bl_bridge') = false;
param_map_group_b('bl_rollover') = true;
fsi_abz_fibs_ds_vecsv_testmain(param_map_group_a)

%% *Case A2*, same as A1, but fl_r_inf (0.10) > fl_r_fbr (0.6)
% Different Formal and informal rate

param_map_group_b('bl_bridge') = false;
param_map_group_b('bl_rollover') = true;
fsi_abz_fibs_ds_vecsv_testmain(param_map_group_b);

%% *Case B1* No Default, Informal Bridge Loan, Formal Menu (formal no rollover), fl_r_inf (0.07) > fl_r_fbr (0.6)
% Same Formal and informal rate

param_map_group_a('bl_bridge') = true;
fsi_abz_fibs_ds_vecsv_testmain(param_map_group_a)

%% *Case B2*, same as B1, but fl_r_inf (0.10) > fl_r_fbr (0.6)
% Different Formal and informal rate

param_map_group_b('bl_bridge') = true;
fsi_abz_fibs_ds_vecsv_testmain(param_map_group_b)
