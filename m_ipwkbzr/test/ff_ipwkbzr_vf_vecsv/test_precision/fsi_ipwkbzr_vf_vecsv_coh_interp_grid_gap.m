%% Tests the IPKWBZ_VF_VECSV Algorithm with varying fl_coh_interp_grid_gap
% For benchmark simulation fl_coh_interp_grid_gap = 0.025;
% Time first it_param_set = 3, then show show results

close all

ar_fl_coh_interp_grid_gap = [0.2, 0.025, 0.01];
% ar_it_w_n = [25, 50];

for fl_coh_interp_grid_gap = ar_fl_coh_interp_grid_gap

    %% Simulate with current ar_it_w
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_coh_interp_grid_gap = ' num2str(fl_coh_interp_grid_gap)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    it_param_set = 4;
    [param_map, support_map] = ffs_ipwkbzr_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_w_perc_n') = 100;
    param_map('it_ak_perc_n') = param_map('it_w_perc_n');    
    param_map('it_z_n') = 15;
    
    param_map('fl_coh_interp_grid_gap') = fl_coh_interp_grid_gap;
    param_map('fl_w_interp_grid_gap') = fl_coh_interp_grid_gap;
    
    param_map('it_c_interp_grid_gap') = 10^-4;
    

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    % Call Program
    ff_ipwkbzr_vf_vecsv(param_map, support_map);


end
