%% Tests the IKWZ_VF_VECSV Algorithm with varying fl_coh_interp_grid_gap
% For benchmark simulation fl_coh_interp_grid_gap = 0.025;
% Time first it_param_set = 3, then show show results

close all

ar_fl_coh_interp_grid_gap = [0.025, 0.0125, 0.00675];
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
    [param_map, support_map] = ffs_akz_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_w_n') = 750;
    param_map('it_ak_n') = param_map('it_w_n');
    param_map('it_z_n') = 11;
    param_map('fl_coh_interp_grid_gap') = fl_coh_interp_grid_gap;
    param_map('it_c_interp_grid_gap') = 10^-4;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    % Call Program
    ff_iwkz_vf_vecsv(param_map, support_map);


end
