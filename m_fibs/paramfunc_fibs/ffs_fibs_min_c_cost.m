%% Formal Informal Cost Minimizer
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [ar_max_c_nobridge, ar_inf_borr_nobridge, ar_for_borr, ar_for_save] = ffs_fibs_min_c_cost(varargin)
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
% @param bl_b_is_principle boolean solving with aggregate savings as
% savings + debt principles + interests, or just principles no interests.
% if true, principels only, no interests.
%
% @param fl_r_inf float informal interest rate
%
% @param fl_r_fsv float (formal) savings interest rate
%
% @param fl_r_fbr float borrowing interest rate
%
% @param ar_aprime array N by 1 level of aggregate borrowing excluding
% bridge loan. Note that bridge loan is needed if coh is negative and
% households can not pay back principle and interests. This must be
% negative.
%
% @return ar_max_c_nobridge array N by 1 next period consumption cost
% (_bl_b_is_principle_ == true), or this period consumption gain
% (_bl_b_is_principle_ == false) based on choosing optimally between formal
% and informal, borrowing and savings joint categories, given either total
% borrowing in principles or principles + interest rate from ar_aprime.
%
% @return ar_inf_borr_nobridge array N by 1 informal borrowing choices
% (Excluding Informal Bridge loans, calculated elsewhere) which could come
% from informal borrowing only if that minimizes consumption cost, or joint
% formal and informal borrowing if that is the cost minimizing choice. if
% _bl_b_is_principle_ == true, then this includes just the principles,  no
% intrest rates. if _bl_b_is_principle_ == false, that means this includes
% interest rates costs as well as principles costs.
%
% @return ar_for_borr array N by 1 formal borrowing choice that minimizes
% consumption costs given fixed _ar_aprime_. Could come from formal
% borrowing alone (which shows up as joint formal and something else where
% the other choice is 0), or formal + informal joint borrow, or formal
% borrowing and formal savings.
%
% @return ar_for_save array N by 1 this is the formal savings choice when
% households are borrowing. Households coulds save just for savings, no
% borrowing as well, that is not captured here.
%
% @example
%
%   bl_input_override = true;
%   [ar_max_c_nobridge, ar_inf_borr_nobridge, ar_for_borr, ar_for_save] = ...
%        ffs_akz_set_functions(bl_b_is_principle, fl_r_inf, fl_r_fsv, ...
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
    [bl_b_is_principle, fl_r_inf, fl_r_fsv, ...
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
    params_group = values(param_map, {'bl_b_is_principle', 'fl_r_inf', 'fl_r_fsv'});
    [bl_b_is_principle, fl_r_inf, fl_r_fsv] = params_group{:};
    
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

ar_max_c_nobridge = zeros(size(ar_aprime_nobridge));
ar_inf_borr_nobridge = zeros(size(ar_aprime_nobridge));
ar_for_borr = zeros(size(ar_aprime_nobridge));
ar_for_save = zeros(size(ar_aprime_nobridge));

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
    ar_for_save(ar_aprime_nobridge_pos_idx) = ar_aprime_nobridge(ar_aprime_nobridge_pos_idx);

    if (bl_b_is_principle)

        % consumption gain next period due to savings
        ar_max_c_nobridge(ar_aprime_nobridge_pos_idx) = ...
            ar_aprime_nobridge(ar_aprime_nobridge_pos_idx)*(1+fl_r_fsv);

    else

        % consumption loss today due to savings
        ar_max_c_nobridge(ar_aprime_nobridge_pos_idx) = ...
            (-1)*ar_aprime_nobridge(ar_aprime_nobridge_pos_idx)./(1+fl_r_fsv);

    end
    
end

%% Proceed to Process Borrowing if aprime Array has borrowing

if (sum(~ar_aprime_nobridge_pos_idx))
    
    %% Compute Consumption Informal Borrowing only
    
    % Generate c Vectors
    if (bl_b_is_principle)
        
        % c_infonly is cost of borrowing in next period consumption
        b_infonly_inf = (ar_aprime_nobridge_br);
        c_infonly = b_infonly_inf.*(1+fl_r_inf);
        
    else
        
        % c_infonly is the gain from borrowing in current period consumption
        b_infonly_inf = (ar_aprime_nobridge_br);
        c_infonly = (-1).*b_infonly_inf./(1+fl_r_inf);
        
    end
    
    % Display
    if (bl_display_minccost)
        tab_c_infonly = table(ar_aprime_nobridge_br, b_infonly_inf, c_infonly);
        disp(['informal borrow interest (fl_r_inf):', num2str(fl_r_inf)]);
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
        b_infforb_inf = (ar_aprime_nobridge_br - ar_a_grid_ceil_principle);
        b_infforb_for = (ar_a_grid_ceil_principle);
        c_infforb =  (b_infforb_inf.*(1+fl_r_inf) + ar_a_grid_ceil_wthr);
        
    else
        
        % c_infforb is the gain from borrowing in current period consumption
        b_infforb_inf = ar_aprime_nobridge_br - ar_a_grid_ceil_wthr;
        b_infforb_for = ar_a_grid_ceil_wthr;
        c_infforb = -1*((b_infforb_inf./(1+fl_r_inf) + ar_a_grid_ceil_principle));
        
    end
    
    % Display
    if (bl_display_minccost)
        tab_c_infforb = table(ar_aprime_nobridge_br, b_infforb_for, b_infforb_inf, c_infforb);
        disp(['formal block size (ar_forbrblk):', num2str(ar_forbrblk)]);
        disp(['formal borrow interest (ar_forbrblk_r):', num2str(ar_forbrblk_r)]);
        disp(['informal borrow interest (fl_r_inf):', num2str(fl_r_inf)]);
        disp(tab_c_infforb);
    end
    
    %% Compute Consumption Values Formal Borrowing + Formal Savings
    
    % Generate c Vectors
    if (bl_b_is_principle)
        
        % c_forbrsv is cost of borrowing in next period consumption
        b_forbrsv_sav = (ar_aprime_nobridge_br - ar_a_grid_floor_principle);
        b_forbrsv_brr = (ar_a_grid_floor_principle);
        c_forbrsv = (b_forbrsv_sav*(1+fl_r_fsv) + ar_a_grid_floor_wthr);
        
    else
        
        % c_forbrsv is the gain from borrowing in current period consumption
        b_forbrsv_sav = ar_aprime_nobridge_br - ar_a_grid_floor_wthr;
        b_forbrsv_brr = ar_a_grid_floor_wthr;
        c_forbrsv = -1*((b_forbrsv_sav./(1+fl_r_fsv) + ar_a_grid_floor_principle));
        
    end
    
    % if b_forbrsv_sav < 0, largest formal borrow grid not large enough. set to
    % 0 so addition of different formal and informal choice categories work.
    b_forbrsv_brr(b_forbrsv_sav < 0) = 0;
    c_forbrsv(b_forbrsv_sav < 0) = nan;
    b_forbrsv_sav(b_forbrsv_sav < 0) = 0;
    
    % Display
    if (bl_display_minccost)
        tab_c_forbrsv = table(ar_aprime_nobridge_br, b_forbrsv_brr, b_forbrsv_sav, c_forbrsv);
        disp(['formal block size (ar_forbrblk):', num2str(ar_forbrblk)]);
        disp(['formal borrow interest (ar_forbrblk_r):', num2str(ar_forbrblk_r)]);
        disp(['savings interest (fl_r_fsv):', num2str(fl_r_fsv)]);
        disp(tab_c_forbrsv);
    end
    
    %% Maximize Consumption For non-Bridge Loan Component
    % Trying to minimize consumption costs
    
    if (bl_b_is_principle)
        
        % minimize consumption cost to the future (also maximizing)
        [ar_max_c_nobridge_br, max_idx] = max([c_infonly c_infforb c_forbrsv],[], 2);
        
    else
        
        % maximum consumption today
        [ar_max_c_nobridge_br, max_idx] = max([c_infonly c_infforb c_forbrsv],[], 2);
        
    end
    
    if (bl_display_minccost)
        tab_c_max = table(ar_aprime_nobridge_br, ar_max_c_nobridge_br, max_idx, c_infonly, c_infforb, c_forbrsv);
        disp(tab_c_max);
    end
    
    %% Borrowing count up borrowing from different sources
    % Informal borrowing comes from informal only or inffor both

    % Consumption when borrowing
    ar_max_c_nobridge(~ar_aprime_nobridge_pos_idx) = ar_max_c_nobridge_br;

    % Informal Borrowing
    ar_inf_borr_nobridge(~ar_aprime_nobridge_pos_idx) = ...
        b_infonly_inf.*(max_idx == 1) + b_infforb_inf.*(max_idx == 2);

    % Formal Borrowing
    ar_for_borr(~ar_aprime_nobridge_pos_idx) = ...
        b_infforb_for.*(max_idx == 2) + b_forbrsv_brr.*(max_idx == 3);

    % Formal Savings
    ar_for_save(~ar_aprime_nobridge_pos_idx) = ...
        b_forbrsv_sav.*(max_idx == 3);
    
end


%% Display Final 
if (bl_display_minccost)
    
    if (bl_b_is_principle)
        ar_average_r = (-1)*(ar_max_c_nobridge./ar_aprime_nobridge);
    else
        ar_average_r = (-1)*(ar_aprime_nobridge./ar_max_c_nobridge);
    end
    
    tab_opti_borrow = table(ar_aprime_nobridge, ar_max_c_nobridge, ...
        ar_average_r, ar_inf_borr_nobridge, ar_for_borr, ar_for_save);
    disp(tab_opti_borrow);
end


end
