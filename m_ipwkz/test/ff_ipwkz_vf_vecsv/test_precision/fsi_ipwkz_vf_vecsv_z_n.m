%% Tests the IPKWZ_VF_VECSV Algorithm with varyin it_z_n
% For benchmark simulation it_z_n = 15
% Time first it_param_set = 3, then show show results

close all

ar_it_z_n = [9, 15, 27];
% ar_it_w_n = [25, 50];

for it_z_n = ar_it_z_n

    %% Simulate with current ar_it_w
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['it_z_n = ' num2str(it_z_n)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    it_param_set = 4;
    [param_map, support_map] = ffs_ipwkz_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_w_perc_n') = 100;
    param_map('it_ak_perc_n') = param_map('it_w_perc_n');    
    param_map('it_z_n') = it_z_n;
    
    param_map('fl_coh_interp_grid_gap') = 0.1;
    param_map('fl_w_interp_grid_gap') = 0.1;    
    param_map('it_c_interp_grid_gap') = 10^-4;
    
    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    % Call Program
    ff_ipwkz_vf_vecsv(param_map, support_map);


end

