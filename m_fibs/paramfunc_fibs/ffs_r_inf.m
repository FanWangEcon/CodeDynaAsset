%% Generate Informal Borrowing Interest Rates
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
function [mt_r_inf_trans, ar_r_inf_stationary, ar_r_inf] = ffs_r_inf(varargin)
%% FFS_R_INF informal interest rate distribution
% Households draw informal interest rates. Probability for each discrete
% element in the sample space of the discrete random variable for informal
% interest rate could be estimated.
%
%

%% Default

bl_input_override = 0;
if (length(varargin) == 3)
    bl_input_override = varargin{3};
end
if (bl_input_override)
    % override when called from outside
    [param_map, support_map, ~] = varargin{:};
else
    % default internal run
    [param_map, support_map] = ffs_fibs_set_default_param();
    default_maps = {param_map, support_map};
    
    bl_display_defparam = true;
    
    % numvarargs is the number of varagin inputted
    [default_maps{1:length(varargin)}] = varargin{:};
    param_map = [param_map; default_maps{1}];
    support_map = [support_map; default_maps{2}];
end

%% Parse Parameters


%% IID

ar_r_inf = [0, 0.05, 0.10, 0.20];
ar_r_inf_stationary = [0.25, 0.25, 0.25, 0.25];
mt_r_inf_trans = repmat(ar_r_inf_stationary, [length(ar_r_inf_stationary), 1]);

if (bl_display_defparam)
    disp('ar_r_inf');
    disp(ar_r_inf);
    disp('ar_r_inf_stationary');
    disp(ar_r_inf_stationary);
    disp('mt_r_inf_trans');
    disp(mt_r_inf_trans);
end

end
