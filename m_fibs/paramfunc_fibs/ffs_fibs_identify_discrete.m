%% Identify Discrete Choices from Continuous Formal Informal Choices
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%
% @seealso
%
% * Formal Borrowing Grid: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_for_br_block_gen.html ffs_for_br_block_gen>
% * Informal Interest Rates: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_r_inf.html ffs_r_inf>
% * Match Borrowing to Formal Grid: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_for_br_block_match.html ffs_for_br_block_match>
% * Optimize Formal and Informal, Borrowing and Savings Joint Choices: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_min_c_cost.html ffs_fibs_min_c_cost>
% * Bridge Loan: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_inf_bridge.html ffs_fibs_inf_bridge>
% * Overall Optimization: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_min_c_cost_bridge.html ffs_fibs_min_c_cost_bridge>
% * Discrete Choices: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_identify_discrete.html ffs_fibs_identify_discrete>
%

%%
function [result_map] = ffs_fibs_identify_discrete(varargin)
%% FFS_FIBS_IDENTIFY_DISCRETE generates categorical formal and informal from cts
% After solving model, have continuous formal and informal choices,
% generate from these categorical outcomes various, append categorical
% outcomes to existing result_map. 
%

%% Default 

bl_input_override = 0;
if (length(varargin) == 2)
    bl_input_override = varargin{2};
end

if (bl_input_override)
    
    % override when called from outside
    [result_map, ~] = varargin{:};
    bl_display_fibs_identify_discrete = false;
    
else
    
    close all
    
    % Initiate
    result_map = containers.Map('KeyType','char', 'ValueType','any');
    
    % Random Matrixes
    rng(123)
    it_n_rows = 5;
    it_n_cols = 5;
    
    rand_zeros_a = round(rand([it_n_rows, it_n_cols]));
    rand_zeros_b = round(rand([it_n_rows, it_n_cols]));
    rand_zeros_c = round(rand([it_n_rows, it_n_cols]));
    rand_zeros_d = round(rand([it_n_rows, it_n_cols]));
    
    mt_pol_b_bridge = -round(rand([it_n_rows, it_n_cols])*20);
    mt_pol_b_bridge(rand_zeros_a == 0) = 0;
    
    mt_pol_inf_borr_nobridge = -round(rand([it_n_rows, it_n_cols])*5);
    mt_pol_inf_borr_nobridge(rand_zeros_b == 0) = 0;

    mt_pol_for_borr = -round(rand([it_n_rows, it_n_cols])*20);
    mt_pol_for_borr(rand_zeros_c == 0) = 0;

    mt_pol_for_save = round(rand([it_n_rows, it_n_cols])*20);
    mt_pol_for_save(rand_zeros_d == 0) = 0;
    
    result_map('cl_mt_pol_b_bridge') = {mt_pol_b_bridge, zeros(1)};
    result_map('cl_mt_pol_inf_borr_nobridge') = {mt_pol_inf_borr_nobridge, zeros(1)};
    result_map('cl_mt_pol_for_borr') = {mt_pol_for_borr, zeros(1)};
    result_map('cl_mt_pol_for_save') = {mt_pol_for_save, zeros(1)};
    
    % Identify
    bl_display_fibs_identify_discrete = true;
end

%% Parse Parameters

params_group = values(result_map, {'cl_mt_pol_b_bridge', 'cl_mt_pol_inf_borr_nobridge', ...
    'cl_mt_pol_for_borr', 'cl_mt_pol_for_save'});
[cl_mt_pol_b_bridge, cl_mt_pol_inf_borr_nobridge, cl_mt_pol_for_borr, cl_mt_pol_for_save] = params_group{:};
[mt_pol_b_bridge, mt_pol_inf_borr_nobridge, mt_pol_for_borr, mt_pol_for_save] = ...
    deal(cl_mt_pol_b_bridge{1}, cl_mt_pol_inf_borr_nobridge{1}, cl_mt_pol_for_borr{1}, cl_mt_pol_for_save{1});

%% Generate Categories

% Generate Binary Outcomes
mt_it_for_borr_idx = (mt_pol_for_borr ~= 0);
mt_it_for_save_idx = (mt_pol_for_save ~= 0);
mt_it_inf_borr_nobridge_idx = (mt_pol_inf_borr_nobridge ~= 0);
mt_it_b_bridge_idx = (mt_pol_b_bridge ~= 0);    

% Generate Multinomial Outcomes
mt_it_for_only_nbdg = ( mt_it_for_borr_idx & ~mt_it_for_save_idx & ~mt_it_inf_borr_nobridge_idx);
mt_it_inf_only_nbdg = (~mt_it_for_borr_idx & ~mt_it_for_save_idx &  mt_it_inf_borr_nobridge_idx);
mt_it_frin_brr_nbdg = ( mt_it_for_borr_idx & ~mt_it_for_save_idx &  mt_it_inf_borr_nobridge_idx);
mt_it_fr_brrsv_nbdg = ( mt_it_for_borr_idx &  mt_it_for_save_idx & ~mt_it_inf_borr_nobridge_idx);
mt_it_frmsavng_only = (~mt_it_for_borr_idx &  mt_it_for_save_idx & ~mt_it_inf_borr_nobridge_idx);

%% Appending Results to result_map

result_map('mt_it_b_bridge_idx') = mt_it_b_bridge_idx;

result_map('mt_it_for_only_nbdg') = mt_it_for_only_nbdg;
result_map('mt_it_inf_only_nbdg') = mt_it_inf_only_nbdg;
result_map('mt_it_frin_brr_nbdg') = mt_it_frin_brr_nbdg;
result_map('mt_it_fr_brrsv_nbdg') = mt_it_fr_brrsv_nbdg;
result_map('mt_it_frmsavng_only') = mt_it_frmsavng_only;

%% Display
if (bl_display_fibs_identify_discrete)

    disp('mt_pol_b_bridge');
    disp(mt_pol_b_bridge);
    disp('mt_pol_inf_borr_nobridge');
    disp(mt_pol_inf_borr_nobridge);
    disp('mt_pol_for_borr');
    disp(mt_pol_for_borr);
    disp('mt_pol_for_save');
    disp(mt_pol_for_save);    
    
    disp('mt_pol_b_bridge');
    disp(mt_pol_b_bridge);
    disp('mt_it_b_bridge_idx');
    disp(mt_it_b_bridge_idx);

    disp('mt_pol_for_borr');
    disp(mt_pol_for_borr);
    disp('mt_it_for_only_nbdg');
    disp(mt_it_for_only_nbdg);
    
    disp('mt_pol_inf_borr_nobridge');
    disp(mt_pol_inf_borr_nobridge);
    disp('mt_it_inf_only_nbdg');
    disp(mt_it_inf_only_nbdg);
    
    disp('mt_pol_for_borr');
    disp(mt_pol_for_borr);
    disp('mt_pol_inf_borr_nobridge');
    disp(mt_pol_inf_borr_nobridge);
    disp('mt_it_frin_brr_nbdg');
    disp(mt_it_frin_brr_nbdg);
    
    disp('mt_pol_for_borr');
    disp(mt_pol_for_borr);
    disp('mt_pol_for_save');
    disp(mt_pol_for_save);
    disp('mt_it_fr_brrsv_nbdg');
    disp(mt_it_fr_brrsv_nbdg);
    
    disp('mt_pol_for_save');
    disp(mt_pol_for_save);
    disp('mt_it_frmsavng_only');
    disp(mt_it_frmsavng_only);
    
    disp('sum of included discrete categories');
    disp(mt_it_for_only_nbdg + mt_it_inf_only_nbdg + mt_it_frin_brr_nbdg + ...
         mt_it_fr_brrsv_nbdg + mt_it_frmsavng_only);
    
end

end
