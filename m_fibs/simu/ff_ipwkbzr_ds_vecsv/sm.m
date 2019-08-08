%% Simulate the Model Given Parameters, Export Model Outcomes to File to be jointly analyzed with data
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
% sm stands for simu match. 

%% Generate Model Simulation Results
clear all;
close all;

it_param_set = 9;
[param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);
support_map('st_mat_path') = "C:\Users\fan\Documents\Dropbox (UH-ECON)\Project_Nguyen_Townsend_Wang_ThaiFinChoice\data_simu\";
support_map('st_mat_prefix') = "ipwkbzrb_";
support_map('st_mat_name_main') = "sm";
support_map('st_mat_suffix') = ["_" num2str(it_param_set)];

[armt_map, func_map] = ffs_ipwkbzr_fibs_get_funcgrid(param_map, support_map);
result_map = ff_ipwkbzr_fibs_vf_vecsv(param_map, support_map, armt_map, func_map);
result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map);

%% Parse Parameters 

params_group = values(support_map, {'st_mat_path', ...
    'st_mat_prefix', 'st_mat_name_main', 'st_mat_suffix'});
[st_mat_path, st_mat_prefix, st_mat_name_main, st_mat_suffix] = params_group{:};

%% States Space Arrays

ar_a = armt_map('ar_a');
ar_z_r_infbr_mesh_wage_w1r2 = armt_map('ar_z_r_infbr_mesh_wage_w1r2');
ar_z_wage_mesh_r_infbr_w1r2 = armt_map('ar_z_wage_mesh_r_infbr_w1r2');
ar_z_r_infbr_prob = armt_map('ar_z_r_infbr_prob');
ar_z_r_infbr = armt_map('ar_z_r_infbr');

tb_outcomes = result_map('tb_outcomes');

%% Overall A Choices

cl_mt_pol_a = result_map('cl_mt_pol_a');
ds_stats_pol_a_map = cl_mt_pol_a{2};
ds_stats_pol_a_map_keys = ds_stats_pol_a_map.keys;

ar_pol_a_choice_prob_byY = ds_stats_pol_a_map('ar_choice_prob_byY');
ar_pol_a_choice_unique_sorted_byY = ds_stats_pol_a_map('ar_choice_unique_sorted_byY');
ar_pol_a_choice_percentiles_byY = ds_stats_pol_a_map('ar_choice_percentiles');
ar_pol_a_fl_percentiles_byY = ds_stats_pol_a_map('ar_fl_percentiles');

%% Overall K Choices

cl_mt_pol_k = result_map('cl_mt_pol_k');
ds_stats_pol_k_map = cl_mt_pol_k{2};
ds_stats_pol_k_map_keys = ds_stats_pol_k_map.keys;

ar_pol_k_choice_prob_byY = ds_stats_pol_k_map('ar_choice_prob_byY');
ar_pol_k_choice_unique_sorted_byY = ds_stats_pol_k_map('ar_choice_unique_sorted_byY');
ar_pol_k_choice_percentiles_byY = ds_stats_pol_k_map('ar_choice_percentiles');
ar_pol_k_fl_percentiles_byY = ds_stats_pol_k_map('ar_fl_percentiles');

%% Overall COH Choices

cl_mt_coh = result_map('cl_mt_coh');
ds_stats_coh_map = cl_mt_coh{2};
ds_stats_coh_map_keys = ds_stats_coh_map.keys;

ar_coh_choice_prob_byY = ds_stats_coh_map('ar_choice_prob_byY');
ar_coh_choice_unique_sorted_byY = ds_stats_coh_map('ar_choice_unique_sorted_byY');
ar_coh_choice_percentiles_byY = ds_stats_coh_map('ar_choice_percentiles');
ar_coh_fl_percentiles_byY = ds_stats_coh_map('ar_fl_percentiles');

%% Export To Mat

mkdir(support_map('st_mat_path'));

clear armt_map func_map param_map support_map result_map

st_file_name = cell2mat([st_mat_prefix st_mat_name_main st_mat_suffix]);
save(strcat(st_mat_path, st_file_name));
