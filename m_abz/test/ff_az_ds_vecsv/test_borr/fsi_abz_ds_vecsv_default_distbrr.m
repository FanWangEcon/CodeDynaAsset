%% With Interest Rate Shocks, What is the Distribution of the Borrowing Interest Rate
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% @seealso
%
% * test speed: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_speed/html/fsi_abz_ds_vecsv_speed.html fsi_abz_ds_vecsv_speed>
% * test joint *RANDOM*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_joint/html/fsi_abz_ds_vecsv_joint_rand.html fsi_abz_ds_vecsv_joint_rand>
%
% * test price no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_price/html/fsi_abz_ds_vecsv_price_nbc_cross.html fsi_abz_ds_vecsv_price_nbc_cross>
% * test price default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_price/html/fsi_abz_ds_vecsv_price_default_cross.html fsi_abz_ds_vecsv_price_default_cross>
%
% * test interest rate no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc.html fsi_abz_ds_vecsv_nbc>
% * test interest rate no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc_cross.html fsi_abz_ds_vecsv_nbc_cross>
% * test interest rate no default *GRID*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_nbc_grid.html fsi_abz_ds_vecsv_nbc_grid>
% * test interest rate default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default.html fsi_abz_ds_vecsv_default>
% * test interest rate default *V(a,z)* Comparison: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default_compaz.html fsi_abz_ds_vecsv_default_compaz>
% * test interest rate default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default_cross.html fsi_abz_ds_vecsv_default_cross>
% * test interest rate default *GRID*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_borr/html/fsi_abz_ds_vecsv_default_grid.html fsi_abz_ds_vecsv_default_grid>
%
% * test shock default (very low cmin): <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_default_lowcmin.html fsi_abz_ds_vecsv_shk_default_lowcmin>
% * test shock no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_nbc.html fsi_abz_ds_vecsv_shk_nbc>
% * test shock no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_nbc_cross.html fsi_abz_ds_vecsv_shk_nbc_cross>
% * test shock default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_default.html fsi_abz_ds_vecsv_shk_default>
% * test shock default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_shock/html/fsi_abz_ds_vecsv_shk_default_cross.html fsi_abz_ds_vecsv_shk_default_cross>
%
% * test preference no default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_nbc.html fsi_abz_ds_vecsv_pref_nbc>
% * test preference no default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_nbc_cross.html fsi_abz_ds_vecsv_pref_nbc_cross>
% * test preference default: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_default.html fsi_abz_ds_vecsv_pref_default>
% * test preference default *CROSS*: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_default_cross.html fsi_abz_ds_vecsv_pref_default_cross>
% * test preference default (very low cmin): <https://fanwangecon.github.io/CodeDynaAsset/m_abz/test/ff_az_ds_vecsv/test_pref/html/fsi_abz_ds_vecsv_pref_default_lowcmin.html fsi_abz_ds_vecsv_pref_default_lowcmin>
%

%% Solve the Model
% Wide interest rate difference between min and max, uniform exogenous
% distribution.

clear all;
close all;

it_param_set = 9;
[param_map, support_map] = ffs_abz_set_default_param(it_param_set);
% param_map('it_a_n') = 50;
% param_map('it_z_wage_n') = 15;
param_map('fl_z_r_borr_max') = 0.50;
param_map('fl_z_r_borr_min') = 0.00;
% param_map('st_z_r_borr_drv_prb_type') = 'unif';
param_map('st_z_r_borr_drv_prb_type') = 'unif';
% param_map('fl_z_r_borr_poiss_mean') = 10;
param_map('fl_z_r_borr_n') = 15;
param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

[armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map);
result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);
result_map = ff_az_ds_vecsv(param_map, support_map, armt_map, func_map, result_map);

%% Obtain The Joint Distribution of Shocks and Savings/Borrowing
% Even though the interest rate is exogenous, borrowing and savings
% decisions are endogenous and are a function of the borrowing rate. The
% realized borrowing interest rate distribution is only for households that
% are borrowing. 

ar_a = armt_map('ar_a');
ar_z_r_borr_mesh_wage = armt_map('ar_z_r_borr_mesh_wage');
ar_z_wage_mesh_r_borr = armt_map('ar_z_wage_mesh_r_borr');
ar_z_r_borr_prob = armt_map('ar_z_r_borr_prob');
ar_z_r_borr = armt_map('ar_z_r_borr');

tb_outcomes = result_map('tb_outcomes');
cl_mt_pol_a = result_map('cl_mt_pol_a');
ds_stats_pol_a_map = cl_mt_pol_a{2};
ds_stats_pol_a_map_keys = ds_stats_pol_a_map.keys;

mt_choice_prob_byYZ = ds_stats_pol_a_map('mt_choice_prob_byYZ');
ar_choice_unique_sorted_byY = ds_stats_pol_a_map('ar_choice_unique_sorted_byY');

% disp(mt_choice_prob_byYZ);
if (norm(sum(mt_choice_prob_byYZ,'all') - 1) <= -1e-12)    
    error('error:sum(mt_choice_prob_byYZ,''all'') ~= 1');
end

%% Marginal Distributional of Realized Borrowing Interest Rate
% Adding up across borrowing levels for each borrowing rate shock

it_z_wage_n = param_map('it_z_wage_n');
fl_z_r_borr_n = param_map('fl_z_r_borr_n');
[it_mt_row_n, it_mt_col_n] = size(mt_choice_prob_byYZ);

mt_choice_prob_byYrZ = zeros(it_mt_row_n, fl_z_r_borr_n);
ar_z_r_borr_dup = zeros(1, fl_z_r_borr_n);
for it_z_r_borr_ctr=1:1:fl_z_r_borr_n
    
    it_start_col = it_z_wage_n*(it_z_r_borr_ctr-1) + 1;
    it_end_col = it_start_col + it_z_wage_n - 1;
    
    ar_z_r_borr_dup(it_z_r_borr_ctr) = ar_z_r_borr_mesh_wage(it_start_col);
    mt_choice_prob_byYrZ(:, it_z_r_borr_ctr) = sum(mt_choice_prob_byYZ(:, it_start_col:it_end_col),2);
end


% R must equal to what was set
if (ar_z_r_borr_dup ~= ar_z_r_borr)
    error('error:ar_z_r_borr_dup ~= ar_z_r_borr');
end

% Sum Up Marginal Probability f(z), must equal to 
if (norm(sum(mt_choice_prob_byYrZ,1)- ar_z_r_borr_prob) <= -1e-12)
    error('error:sum(mt_choice_prob_byYrZ,1) ~= ar_z_r_borr_prob');
end

if (norm(sum(mt_choice_prob_byYrZ,'all') - 1) <= -1e-12)
    error('error:sum(mt_choice_prob_byYrZ,''all'') ~= 1');
end

%% Frac of Borrower for each R: P(a <= 0 | z)

ar_it_choice_unique_sorted_byY_neg = (ar_choice_unique_sorted_byY < 0);
mt_choice_prob_byYrZ_borr = mt_choice_prob_byYrZ(ar_it_choice_unique_sorted_byY_neg,:);

% ar_p_a_le_zr_r_borr = sum(P(a, z) 1(a <=0))
ar_p_a_le_zr_r_borr = sum(mt_choice_prob_byYrZ_borr, 1);

% ar_p_a_le_zr_r_borr = sum(P( z | a <=0))
ar_p_r_borr = ar_p_a_le_zr_r_borr/sum(ar_p_a_le_zr_r_borr);

% ar_e_a_le_zr_r_borr = sum(a*P(a, z) 1(a <=0))
ar_Ebr_le_zr_r_borr = sum((ar_choice_unique_sorted_byY(ar_choice_unique_sorted_byY < 0) ...
                          + zeros(1,fl_z_r_borr_n)) .* mt_choice_prob_byYrZ_borr, 1);
                      
% share of total borrowing by r
ar_Ebr_share_r_borr = ar_Ebr_le_zr_r_borr/sum(ar_Ebr_le_zr_r_borr);

                      
disp(ar_z_r_borr_dup);
disp(ar_z_r_borr_prob);
disp(ar_p_r_borr);
disp(ar_Ebr_share_r_borr);

%% Graph 

figure();
hold on;
bar(ar_z_r_borr', [ar_z_r_borr_prob;ar_p_r_borr;ar_Ebr_share_r_borr]');
legend({'exogenous borrow rate distribution', 'borrow rate distribution conditional on borrowing', 'borrow volumn conditional on rate'}, ...
        'Location', 'northeast', 'color', 'none');
title({'Model With Exogenous Uniform Interest Rate Shocks', 'Borrowing is Endogenous, Realized R dist. differs from exo rate dist'})
ylabel('Fractions')
xlabel('Borrowing Interest Rates')
grid on;
grid minor;