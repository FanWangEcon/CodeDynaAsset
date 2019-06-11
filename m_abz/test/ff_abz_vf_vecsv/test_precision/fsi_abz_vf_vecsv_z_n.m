%% Tests the ABZ_VF_VECSV Algorithm with varyin it_z_n
% For benchmark simulation it_z_n = 15

close all

ar_it_z_n = [9, 15, 27];

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
    [param_map, support_map] = ffs_abz_set_default_param(it_param_set);

    % Simulation Accuracy
    param_map('it_a_n') = 750;
    param_map('it_z_n') = it_z_n;

    % Display Parameters
    support_map('bl_display') = false;
    support_map('bl_display_final') = false;
    support_map('bl_time') = true;
    support_map('bl_profile') = false;

    % Call Program
    ff_abz_vf_vecsv(param_map, support_map);


end
