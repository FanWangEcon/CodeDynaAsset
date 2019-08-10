%% Single Simulate Dense
% Single Dense Simulation to See policy functions clearl 

%% Simulate Set Parameters Defalt
clear all;
close all;

it_param_set = 8;
[param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);

%% Simulate Set Parameters Dense
st_param_which = 'dense';

if (ismember(st_param_which, ["default"]))
    
%     support_map('it_display_every') = 10;
    param_map('it_maxiter_val') = 40;
    
elseif (ismember(st_param_which, ["dense"]))
    support_map('it_display_every') = 1;
    
    param_map('it_maxiter_val') = 40;
    param_map('it_w_perc_n') = 200;
    param_map('it_ak_perc_n') = param_map('it_w_perc_n');
    param_map('it_coh_bridge_perc_n') = param_map('it_w_perc_n')/3;

    param_map('fl_coh_interp_grid_gap') = 0.02;
    param_map('it_c_interp_grid_gap') = 10^-4;
    param_map('fl_w_interp_grid_gap') = 0.02;

    param_map('it_z_wage_n') = 7;
    param_map('fl_z_r_infbr_n') = 5;
    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');
    
end

%% Simulate
ff_ipwkbzr_fibs_ds_wrapper(param_map, support_map);