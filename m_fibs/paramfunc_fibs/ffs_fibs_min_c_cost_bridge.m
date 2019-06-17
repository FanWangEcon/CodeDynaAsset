%% Formal and Informal Choices Given Aprime and cash-on-hand
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [fl_max_c, fl_b_bridge, fl_inf_borr_nobridge, fl_for_borr, fl_for_save] = ...
            ffs_fibs_min_c_cost_bridge(varargin)
%% FFS_FIBS_MIN_C_COST_INF_BRIDGE combine ffs_fibs_min_c_cost + inf_bridge
% Given coh and aprime choice, what are the optimal formal and informal
% joint choices including bridge choices. This function is invoked after
% optimal a-prime choices have been found. This is invoked so that during
% solution, do not have to store these choices. This is the same material
% as what is in ff_abz_fibs_vf.m
%
% @param fl_ap float aprime choice, for example optimal aprime solved
%
% @param fl_coh float cash-on-hand for the aprime choice
%
% @param ar_aprime array N by 1 level of aggregate borrowing excluding
% bridge loan. Note that bridge loan is needed if coh is negative and
% households can not pay back principle and interests. This must be
% negative.
%
% @return fl_max_c float next period consumption cost
% (_bl_b_is_principle_ == true), or this period consumption gain
% (_bl_b_is_principle_ == false) based on choosing optimally between formal
% and informal, borrowing and savings joint categories, given either total
% borrowing in principles or principles + interest rate from ar_aprime.
%
% @return fl_b_bridge float bridge loan debt to pay for unpaid uncovered
% cash-on-hand
%
% @return fl_inf_borr_nobridge float informal borrowing choices
% (Excluding Informal Bridge loans, calculated elsewhere) which could come
% from informal borrowing only if that minimizes consumption cost, or joint
% formal and informal borrowing if that is the cost minimizing choice. if
% _bl_b_is_principle_ == true, then this includes just the principles,  no
% intrest rates. if _bl_b_is_principle_ == false, that means this includes
% interest rates costs as well as principles costs.
%
% @return fl_for_borr float formal borrowing choice that minimizes
% consumption costs given fixed _ar_aprime_. Could come from formal
% borrowing alone (which shows up as joint formal and something else where
% the other choice is 0), or formal + informal joint borrow, or formal
% borrowing and formal savings.
%
% @return fl_for_save float this is the formal savings choice when
% households are borrowing. Households coulds save just for savings, no
% borrowing as well, that is not captured here.
%
% @example
%
%   bl_input_override = true;
%   [fl_max_c, fl_b_bridge, fl_inf_borr_nobridge, fl_for_borr, fl_for_save] = ...
%        ffs_fibs_min_c_cost_bridge(fl_ap, fl_coh, ...
%            param_map, support_map, armt_map, func_map, bl_input_override);
%


%% Default 

bl_input_override = 0;
if (length(varargin) == 7)
    bl_input_override = varargin{7};
end
if (bl_input_override)
    % override when called from outside
    [fl_ap, fl_coh, param_map, support_map, armt_map, func_map, ~] = varargin{:};
    
else
    
    close all
    
    % Default
    it_param_set = 4;
    bl_input_override = true;
    [param_map, support_map] = ffs_abz_fibs_set_default_param(it_param_set);
    [armt_map, func_map] = ffs_abz_fibs_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
    
    % Testing COH and Aprime Vectors
    fl_ap = -10;
    fl_coh = 5;
        
    % Set Display Control
    support_map('bl_display_infbridge') = true;
    support_map('bl_display_minccost') = true;
    
end

%% Parse Parameters

% Gather Inputs from armt_map
params_group = values(armt_map, {'ar_forbrblk', 'ar_forbrblk_r'});
[ar_forbrblk, ar_forbrblk_r] = params_group{:};

% Gather Inputs from param_map
params_group = values(param_map, {'bl_default', 'bl_b_is_principle', 'fl_r_inf', 'fl_r_fsv', 'fl_c_min'});
[bl_default, bl_b_is_principle, fl_r_inf, fl_r_fsv, fl_c_min] = params_group{:};

% func_map
params_group = values(func_map, {'f_cons_coh_fbis', 'f_cons_coh_save'});
[f_cons_coh_fbis, f_cons_coh_save] = params_group{:};

% support_map
params_group = values(support_map, {'bl_display_minccost', 'bl_display_infbridge'});
[bl_display_minccost, bl_display_infbridge] = params_group{:};

%% Compute Consumption given Borrowing and Savings
% find the today's consumption maximizing formal and informal
% choices given a' and coh. The formal and informal choices need to
% generate exactly a', but depending on which formal and informal
% joint choice is used, the consumption cost today a' is different.
% Note here, a is principle + interests. Three areas:
%
% * *CASE A* a' > 0: savings, do not need to optimize over formal and
% informal choices
% * *CASE B* a' < 0 & coh < 0: need bridge loan to pay for unpaid debt, and
% borrowing over-all, need to first pick bridge loan to pay for
% debt, if bridge loan is insufficient, go into default. After
% bridge loan, optimize over formal+informal, borrow+save joint
% choices.
% * *CASE C* $ a' < 0 & coh > 0: do not need to get informal bridge loans,
% optimize over for+inf save, for+save+borr, inf+borr only, for
% borrow only.
%

if (fl_ap < 0)
    
    % Calculate Bridge Loan Borrowing
    if (fl_coh < 0)
        bl_input_override = true;
        [fl_aprime_nobridge, fl_b_bridge, fl_c_bridge] = ffs_fibs_inf_bridge(...
            bl_b_is_principle, fl_r_inf, fl_ap, fl_coh, ...
            bl_display_infbridge, bl_input_override);
        
    else
        
        fl_aprime_nobridge = fl_ap;
        fl_b_bridge = 0;
        fl_c_bridge = 0;
        
    end
    
    % Find Optimal Formal Informal Borrow Save Combo
    % calculate consumption gain from formal + informal
    % borrowing and savings choices.
    bl_input_override = true;
    [fl_max_c_nobridge, fl_inf_borr_nobridge, fl_for_borr, fl_for_save] = ...
        ffs_fibs_min_c_cost(...
        bl_b_is_principle, fl_r_inf, fl_r_fsv, ...
        ar_forbrblk_r, ar_forbrblk, ...
        fl_aprime_nobridge, bl_display_minccost, bl_input_override);
    
    % Compute Consumption given Formal and Informal joint
    % consumption with formal borrow menu + bridge loans.
    fl_max_c = f_cons_coh_fbis(fl_coh, fl_max_c_nobridge + fl_c_bridge);
    
else
    
    % consumption with savings
    fl_max_c = f_cons_coh_save(fl_coh, fl_ap);
    
    % assign values for formal and informal choices
    % possible that fl_coh < 0, but if then fl_ap > 0 is
    % not valid choice
    [fl_b_bridge, fl_inf_borr_nobridge, fl_for_borr, fl_for_save] = deal(0, 0, 0, fl_ap);
    
end

%% Compute Utility With Default
% assign u(c)
if (fl_max_c <= fl_c_min)
    if (bl_default)
        % Replace Consumption if default cmin
        fl_max_c = fl_c_min;
    else
        % Replace Consumption if no default nan
        fl_max_c = 0;
    end
end

end


