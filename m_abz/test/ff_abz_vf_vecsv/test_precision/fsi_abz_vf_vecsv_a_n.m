%% Tests the ABZ_VF_VECSV Algorithm with varyin it_a_n

close all

ar_it_a_n = [100, 750, 2250];
% ar_it_w_n = [25, 50];

for it_a_n = ar_it_a_n

    %% Simulate with current ar_it_w
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['it_a_n = ' num2str(it_a_n)]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('');
    disp('');
    disp('');
    disp('');

    it_param_set = 4;
    [param_map, support_map] = ffs_abz_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_a_n') = it_a_n;
    param_map('it_z_n') = 15;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    % Call Program
    ff_abz_vf_vecsv(param_map, support_map);

end
