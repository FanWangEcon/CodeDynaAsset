%% Test Production Function Parameters (CROSS)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%% Set Shared Parameters

close all;
clear all;

bl_default = true;
ar_fl_z_wage_rho = [0.0, 0.50, 0.99];
ar_fl_z_wage_sig = [0.05, 0.10, 0.3];

ar_fl_z_wage_rho = [0.91464, 0.985];
ar_fl_z_wage_rho = [0.985];
ar_fl_z_wage_sig = [0.20];

%% Simulate Model with schok persistence = 0.0, IID

for fl_z_wage_rho = ar_fl_z_wage_rho
    for fl_z_wage_sig = ar_fl_z_wage_sig

        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp(['fl_z_wage_rho = ' num2str(fl_z_wage_rho)]);
        disp(['fl_z_wage_sig = ' num2str(fl_z_wage_sig)]);
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp('');
        disp('');
        disp('');
        disp('');

        tic;
        bl_input_override = true;
        it_param_set = 9;
        [param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);

        param_map('fl_z_wage_rho') = fl_z_wage_rho;
        param_map('fl_z_wage_sig') = fl_z_wage_sig;
        param_map('bl_default') = bl_default;

        support_map('bl_time') = true;
        support_map('bl_display') = true;
        support_map('it_display_every') = 20;
        support_map('bl_display_defparam') = true;
        support_map('bl_display_final') = true;
        
        support_map('bl_display_dist') = true;        
        support_map('bl_post') = true;
        support_map('bl_display_final_dist') = true;
        
        [armt_map, func_map] = ffs_ipwkbzr_fibs_get_funcgrid(param_map, support_map, bl_input_override);
        result_map = ff_ipwkbzr_fibs_vf_vecsv(param_map, support_map, armt_map, func_map);
        result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
        toc;

        % Snap
        snapnow;

    end
end
