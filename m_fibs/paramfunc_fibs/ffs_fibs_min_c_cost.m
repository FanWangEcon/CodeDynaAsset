%% Formal Informal Cost Minimizer
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [mt_max_c_nobridge, mt_inf_borr_nobridge, mt_for_borr, mt_for_save] = ffs_fibs_min_c_cost(varargin)
%% FFS_FIBS_MIN_C_COST Minimizes to consumption given aprime.
% Suppose I want to borrow 5 dollar today. That will increase current
% consumption by 5 dollars, but its cost to next period cash-on-hand is a
% function of how I get the 5 dollar of borrowing, whether I borrow from
% formal sources, informal or joint sources. Suppose the formal borrowing
% grid has three points, you can borrow 0, 3 or 6, these are the options:
%
% # borrow all 5 informally
% # borrow 3 formally and 2 informally
% # borrow 6 formally and save 1
% # additionally, if 5 is actually on the grid, you would borrow 5 only
% formally if that offers the lower interest rate.
%
% This function here finds which one of of these options is optimal
% conditional on aprime choice. Households choose among the different
% aprime choices given these within aprime choice optimal formal/informal
% allocations.
%
% _ar_aprime_ is the borrowing choices that do not include bridge
% loan. The additional borrowing on top of bridge loans.
%
% The function outputs matrixes in _N_ by _M_ dimensions. _N_ is the
% dimension of the _ar_aprime_ array, and _M_ is the dimension fo the
% _ar_r_inf_ interest rate array. Solves the optimal formal and informal
% joint optimization prblem for *N* by *M* combinations of aggregate
% borrowing levels and informal borrowing interest rate. The idea is that
% there are different aggregate borrowing levels households are choosing
% from, and they face different informal borrowing interest rate shocks. 
%
% @param bl_b_is_principle boolean solving with aggregate savings as
% savings + debt principles + interests, or just principles no interests.
% if true, principels only, no interests.
%
% @param ar_r_inf array *1 by M* array of informal interest rate
%
% @param fl_r_fsv float (formal) savings interest rate
%
% @param fl_r_fbr float borrowing interest rate
%
% @param ar_aprime array *N by 1* level of aggregate borrowing excluding
% bridge loan. Note that bridge loan is needed if coh is negative and
% households can not pay back principle and interests. This must be
% negative.
%
% @return mt_max_c_nobridge array N by M next period consumption cost
% (_bl_b_is_principle_ == true), or this period consumption gain
% (_bl_b_is_principle_ == false) based on choosing optimally between formal
% and informal, borrowing and savings joint categories. This considers both
% interests as well as principles.
%
% @return mt_inf_borr_nobridge array N by M informal borrowing choices
% (Excluding Informal Bridge loans, calculated elsewhere) which could come
% from informal borrowing only if that minimizes consumption cost, or joint
% formal and informal borrowing if that is the cost minimizing choice. ZIf
% _bl_b_is_principle_ == true, then this includes just the principles,  no
% intrest rates. if _bl_b_is_principle_ == false, that means this includes
% interest rates costs as well as principles costs.
%
% @return mt_for_borr array N by M formal borrowing choice that minimizes
% consumption costs given fixed _ar_aprime_. Could come from formal
% borrowing alone (which shows up as joint formal and something else where
% the other choice is 0), or formal + informal joint borrow, or formal
% borrowing and formal savings. If _bl_b_is_principle_ == true, then this
% includes just the principles,  no intrest rates. if _bl_b_is_principle_
% == false, that means this includes interest rates costs as well as
% principles costs.
%
% @return mt_for_save array N by M this is the formal savings choice when
% households are borrowing. Households coulds save just for savings, no
% borrowing as well, that is not captured here. If _bl_b_is_principle_ ==
% true, then this includes just the principles,  no intrest rates. if
% _bl_b_is_principle_ == false, that means this includes interest rates
% costs as well as principles costs.
%
% @example
%
%   bl_input_override = true;
%   [ar_max_c_nobridge, ar_inf_borr_nobridge, ar_for_borr, ar_for_save] = ...
%        ffs_fibs_min_c_cost(bl_b_is_principle, fl_r_inf, fl_r_fsv, ...
%                              ar_forbrblk_r, ar_forbrblk, ...
%                              fl_ap, bl_display_minccost, bl_input_override);
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc/html/ffs_abz_fibs_set_default_param.html ffs_abz_fibs_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc/html/ffs_abz_fibs_get_funcgrid.html ffs_abz_fibs_get_funcgrid>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_for_br_block_match.html ffs_for_br_block_match>
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

%% Default and Parse Parameters

bl_input_override = 0;
if (length(varargin) == 8)
    bl_input_override = varargin{8};
end
if (bl_input_override)
    % override when called from outside
    [bl_b_is_principle, ar_r_inf, fl_r_fsv, ...
        ar_forbrblk_r, ar_forbrblk, ar_aprime_nobridge, bl_display_minccost, ~] = varargin{:};
    
else
    
    close all
    
    % Default
    it_param_set = 4;
    bl_input_override = true;
    [param_map, support_map] = ffs_ipwkbz_fibs_set_default_param(it_param_set);
    
    % Gather Inputs from armt_map
    params_group = values(param_map, ...
        {'fl_r_fbr', 'st_forbrblk_type', 'fl_forbrblk_brmost', 'fl_forbrblk_brleast', 'fl_forbrblk_gap'});
    [fl_r_fbr, st_forbrblk_type, fl_forbrblk_brmost, fl_forbrblk_brleast, fl_forbrblk_gap] = params_group{:};
    [ar_forbrblk, ar_forbrblk_r] = ...
        ffs_for_br_block_gen(fl_r_fbr, st_forbrblk_type, fl_forbrblk_brmost, fl_forbrblk_brleast, fl_forbrblk_gap);
    
    % Gather Inputs from param_map
    params_group = values(param_map, {'bl_b_is_principle', 'fl_r_fsv'});
    [bl_b_is_principle,  fl_r_fsv] = params_group{:};
    
    % Informal Interest Rate Array, works with 1 by 1 or 1 by N
    ar_r_inf = [0.00, 0.085, 0.90];
%     ar_r_inf = [0.085];
    
    % Setting interest rate, if r_inf is very high, that means informal
    % option would generally never be chosen, or options that involve
    % informal options never choice, either here, or outside.
    %     fl_r_inf = 10000;
    
    % Testing COH and Aprime Vectors
    ar_aprime_nobridge = [-20, -14, -11, -6.8, ...
        -5.5, -4.5, -4.1, -1.1, ...
        -0.1, ...
        0.1, 1, 2, 10]';
    
    % Set Display Control
    bl_display_minccost = true;
    
end

%% Initialize Output Arrays

it_N = length(ar_aprime_nobridge);
it_M = length(ar_r_inf);

mt_max_c_nobridge = zeros(it_N, it_M);
mt_inf_borr_nobridge = zeros(it_N, it_M);
mt_for_borr = zeros(it_N, it_M);
mt_for_save = zeros(it_N, it_M);

%% Find if ar_aprime contains Saving
% only one savings option, function here meant for borrowing, but can deal
% with savings as well if they appear.

ar_aprime_nobridge_pos_idx = (ar_aprime_nobridge >= 0);
ar_aprime_nobridge_br = ar_aprime_nobridge(~ar_aprime_nobridge_pos_idx);

%% Get consumption and savings choices if savings
% When households overall save, they could still have had to first pay
% informal lender for bridge loan. For households with positive
% cash-on-hand, if they save, they are only savings, and not borrowing at
% the same time. In practice, some households coule borrow from one sector
% and lend/save in another. This is not a mechanism in this model, but is
% discussed in <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3316939 Wang (2019)>
%

if (sum(ar_aprime_nobridge_pos_idx))
    
    % savings = savings (including or not interest rates)
    mt_for_save(ar_aprime_nobridge_pos_idx) = ar_aprime_nobridge(ar_aprime_nobridge_pos_idx);
    
    if (bl_b_is_principle)
        
        % consumption gain next period due to savings
        mt_max_c_nobridge(ar_aprime_nobridge_pos_idx, :) = ...
            zeros(1,it_M) + ar_aprime_nobridge(ar_aprime_nobridge_pos_idx)*(1+fl_r_fsv);
        
    else
        
        % consumption loss today due to savings
        mt_max_c_nobridge(ar_aprime_nobridge_pos_idx, :) = ...
            zeros(1,it_M) + (-1)*ar_aprime_nobridge(ar_aprime_nobridge_pos_idx)./(1+fl_r_fsv);
        
    end
    
end

%% Proceed to Process Borrowing if aprime Array has borrowing

if (sum(~ar_aprime_nobridge_pos_idx))
    
    %% Compute Consumption Informal Borrowing only
    
    % Generate c Vectors
    if (bl_b_is_principle)
        
        % c_infonly is cost of borrowing in next period consumption
        mt_b_infonly_inf = (ar_aprime_nobridge_br);
        mt_c_infonly = mt_b_infonly_inf.*(1+ar_r_inf);
        
    else
        
        % c_infonly is the gain from borrowing in current period consumption
        mt_b_infonly_inf = (ar_aprime_nobridge_br);
        mt_c_infonly = (-1).*mt_b_infonly_inf./(1+ar_r_inf);
        
    end
    
    % Display
    if (bl_display_minccost)
        tab_c_infonly = table(ar_aprime_nobridge_br, mt_b_infonly_inf, mt_c_infonly);
        disp(['informal borrow interest (fl_r_inf):', num2str(ar_r_inf)]);
        disp(tab_c_infonly);
    end
    
    %% Compute Formal Block Sizes
    % Divide each asset choice by each element of the formal choice grid. The
    % numbers closest to 1 above and below indicates the floor and ceil formal
    % borrowing quantity for joint credit market participation choices.
    
    [ar_a_grid_ceil_principle, ar_a_grid_ceil_wthr, ar_a_grid_floor_principle, ar_a_grid_floor_wthr] = ...
        ffs_for_br_block_match(ar_aprime_nobridge_br, ar_forbrblk, ar_forbrblk_r, bl_b_is_principle);
    
    %% Compute Consumption Values Formal + Informal Borrowing Jointly
    % if right on grid, informal could be zero.
    
    % Generate c Vectors
    if (bl_b_is_principle)
        
        % c_infforb is cost of borrowing in next period consumption
        ar_b_infforb_inf = (ar_aprime_nobridge_br - ar_a_grid_ceil_principle);
        ar_b_infforb_for = (ar_a_grid_ceil_principle);
        mt_c_infforb =  (ar_b_infforb_inf.*(1+ar_r_inf) + ar_a_grid_ceil_wthr);
        
    else
        
        % c_infforb is the gain from borrowing in current period consumption
        ar_b_infforb_inf = ar_aprime_nobridge_br - ar_a_grid_ceil_wthr;
        ar_b_infforb_for = ar_a_grid_ceil_wthr;
        mt_c_infforb = -1*((ar_b_infforb_inf./(1+ar_r_inf) + ar_a_grid_ceil_principle));
        
    end
    
    % Display
    if (bl_display_minccost)
        tab_c_infforb = table(ar_aprime_nobridge_br, ar_b_infforb_for, ar_b_infforb_inf, mt_c_infforb);
        disp(['formal block size (ar_forbrblk):', num2str(ar_forbrblk)]);
        disp(['formal borrow interest (ar_forbrblk_r):', num2str(ar_forbrblk_r)]);
        disp(['informal borrow interest (fl_r_inf):', num2str(ar_r_inf)]);
        disp(tab_c_infforb);
    end
    
    %% Compute Consumption Values Formal Borrowing + Formal Savings
    
    % Generate c Vectors
    if (bl_b_is_principle)
        
        % c_forbrsv is cost of borrowing in next period consumption
        ar_b_forbrsv_sav = (ar_aprime_nobridge_br - ar_a_grid_floor_principle);
        ar_b_forbrsv_brr = (ar_a_grid_floor_principle);
        ar_c_forbrsv = (ar_b_forbrsv_sav*(1+fl_r_fsv) + ar_a_grid_floor_wthr);
        
    else
        
        % c_forbrsv is the gain from borrowing in current period consumption
        ar_b_forbrsv_sav = ar_aprime_nobridge_br - ar_a_grid_floor_wthr;
        ar_b_forbrsv_brr = ar_a_grid_floor_wthr;
        ar_c_forbrsv = -1*((ar_b_forbrsv_sav./(1+fl_r_fsv) + ar_a_grid_floor_principle));
        
    end
    
    % if b_forbrsv_sav < 0, largest formal borrow grid not large enough. set to
    % 0 so addition of different formal and informal choice categories work.
    ar_b_forbrsv_brr(ar_b_forbrsv_sav < 0) = 0;
    ar_c_forbrsv(ar_b_forbrsv_sav < 0) = nan;
    ar_b_forbrsv_sav(ar_b_forbrsv_sav < 0) = 0;
    
    % Display
    if (bl_display_minccost)
        tab_c_forbrsv = table(ar_aprime_nobridge_br, ar_b_forbrsv_brr, ar_b_forbrsv_sav, ar_c_forbrsv);
        disp(['formal block size (ar_forbrblk):', num2str(ar_forbrblk)]);
        disp(['formal borrow interest (ar_forbrblk_r):', num2str(ar_forbrblk_r)]);
        disp(['savings interest (fl_r_fsv):', num2str(fl_r_fsv)]);
        disp(tab_c_forbrsv);
    end
    
    %% Maximize Consumption For non-Bridge Loan Component
    
    mt_c_forbrsv = repmat(ar_c_forbrsv, [1, length(ar_r_inf)]);
    [ar_max_c_nobridge_br, ar_max_idx] = max([mt_c_infonly(:) mt_c_infforb(:) mt_c_forbrsv(:)],[], 2);
    mt_max_c_nobridge_br = reshape(ar_max_c_nobridge_br, size(mt_c_infonly));
    mt_max_idx = reshape(ar_max_idx, size(mt_c_infonly));
    
    if (bl_display_minccost)
        disp('ar_r_inf')
        disp(ar_r_inf)
        
        disp('ar_aprime_nobridge''')
        disp(ar_aprime_nobridge')
        
        tab_c_max = table(ar_aprime_nobridge_br, mt_max_c_nobridge_br, mt_max_idx, mt_c_infonly, mt_c_infforb, ar_c_forbrsv);
        disp(tab_c_max);
    end
    
    %% Borrowing count up borrowing from different sources
    % Informal borrowing comes from informal only or inffor both
    
    % Consumption when borrowing
    mt_max_c_nobridge(~ar_aprime_nobridge_pos_idx, :) = mt_max_c_nobridge_br;
    
    % Informal Borrowing
    mt_inf_borr_nobridge(~ar_aprime_nobridge_pos_idx, :) = ...
        mt_b_infonly_inf.*(mt_max_idx == 1) + ar_b_infforb_inf.*(mt_max_idx == 2);
    
    % Formal Borrowing
    mt_for_borr(~ar_aprime_nobridge_pos_idx, :) = ...
        ar_b_infforb_for.*(mt_max_idx == 2) + ar_b_forbrsv_brr.*(mt_max_idx == 3);
    
    % Formal Savings
    mt_for_save(~ar_aprime_nobridge_pos_idx, :) = ...
        ar_b_forbrsv_sav.*(mt_max_idx == 3);
    
end


%% Display Final
if (bl_display_minccost)
    
    if (bl_b_is_principle)
        ar_average_r = (-1)*(mt_max_c_nobridge./ar_aprime_nobridge);
    else
        ar_average_r = (-1)*(ar_aprime_nobridge./mt_max_c_nobridge);
    end
    
    tab_opti_borrow = table(ar_aprime_nobridge, mt_max_c_nobridge, ...
        ar_average_r, mt_inf_borr_nobridge, mt_for_borr, mt_for_save);
    disp(tab_opti_borrow);
end


end
