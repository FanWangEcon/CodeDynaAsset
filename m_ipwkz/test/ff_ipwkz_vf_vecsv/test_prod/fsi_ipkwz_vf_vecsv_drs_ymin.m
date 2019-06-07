%% Tests the IPKWZ_VF_VECSV Algorithm with varyin w (fl_y_min)

close all

ar_fl_y_min = [0, 0.05, 0.15];
% ar_it_w_n = [25, 50];

for fl_y_min = ar_fl_y_min

    %% Simulate with current ar_it_w
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['fl_y_min = ' num2str(fl_y_min)]);
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
    param_map('it_z_n') = 11;

    param_map('fl_coh_interp_grid_gap') = 0.0125;
%     param_map('fl_w_interp_grid_gap') = 0.0125;
    param_map('fl_w_interp_grid_gap') = (param_map('fl_w_max') - param_map('fl_w_min'))/param_map('it_w_perc_n');
    param_map('it_c_interp_grid_gap') = 10^-4;

    % Iteration Limits
    % param_map('it_maxiter_val') = 2;

    % Turn on 2nd stage graphs
    support_map('bl_graph_evf') = false;
    support_map('bl_display_evf') = false;

    % Production Function Parameters
    param_map('fl_w') = fl_y_min;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = true;
    support_map('bl_time') = true;
    % support_map('bl_profile') = false;


    % Call Program
    ff_ipwkz_vf_vecsv(param_map, support_map);


end
