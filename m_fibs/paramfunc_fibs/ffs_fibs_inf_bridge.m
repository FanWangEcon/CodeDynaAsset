%% Formal Informal Bridge Loan Handling
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%%

function [ar_aprime_nobridge, ar_b_bridge, ar_c_bridge] = ffs_fibs_inf_bridge(varargin)
%% FFS_FIBS_INF_BRIDGE Amount of Informal Borrowing Needed as Bridge Loans
% Bridge loan needed to pay for debt that is still unpaid due to
% insufficient cash on hand. Only informal lender or some other lender
% willing to extend this loan offers it
%
% @param bl_b_is_principle boolean solving with aggregate savings as
% savings + debt principles + interests, or just principles no interests.
% if true, principels only, no interests.
%
% @param fl_r_bridge float interest rate for bridge loan
%
% @param fl_r_fsv float (formal) savings interest rate
%
% @param fl_r_fbr float borrowing interest rate
%
% @param ar_aprime array N by 1 level of aggregate borrowing excluding
% bridge loan. Note that bridge loan is needed if coh is negative and
% households can not pay back principle and interests.
%
% @param ar_coh_today array N by 1 the level of cash-on-hand today, when the
% borrowing and savings decisions are made. If this is positive, then
% households freely borrow, do not need bridge loans. If this is negative
% households need to first borrow to meet bridge loan needs.
%
% @return ar_aprime_nobridge array next period asset choice without debt
% incurred for bridge loans.
%
% @return ar_b_bridge array bridge loan debt to pay for unpaid uncovered
% cash-on-hand
%
% @return ar_c_bridge array consumption gain from bridge loan, or
% consumption costs of bridge loans
%
% @example
%
%   bl_input_override = true;
%   [ar_aprime_nobridge, ar_c_bridge] = ...
%        ffs_fibs_inf_bridge(bl_b_is_principle, fl_r_bridge, ...
%                              ar_aprime, ar_coh_today, ...
%                              bl_display_infbridge, bl_input_override);
%


%% Default and Parse Parameters

bl_input_override = 0;
if (length(varargin) == 6)
    bl_input_override = varargin{6};
end

if (bl_input_override)
    % override when called from outside
    [bl_b_is_principle, fl_r_bridge, ar_aprime, ar_coh_today, bl_display_infbridge, ~] = varargin{:};
else
    close all

    % Default
    it_param_set = 4;
    [param_map, ~] = ffs_abz_fibs_set_default_param(it_param_set);

    % Gather Inputs from param_map
    params_group = values(param_map, {'bl_b_is_principle', 'fl_r_inf'});
    [bl_b_is_principle, fl_r_inf] = params_group{:};

    % For benchmark, assume that the informal lender
    fl_r_bridge = fl_r_inf;

    % Testing COH and Aprime Vectors
    ar_aprime =    [-20. -5,-5, -4.5,-4.5, -0.1,-0.1]';
    ar_coh_today = [-19, 1, -1,   1,-1,      1,0 ]';

    % Set Display Control
    bl_display_infbridge = true;
end

%% Bridge Loan Required
% when coh <= cmin/0. Income does not fully repay debts. Suppose formal
% lenders have strict rules about not allowing for roll-over. Then when
% this happens, households would go to default state if default is allowed.
% If default is not allowed, households would never borrow such that coh
% ends up lower than 0. But now informal lender comes in and is willing to
% offer bridge loan. This bridge loan might be a fraction of total debt
% taken out last period, and it will become share of the debt carried on
% today. Or households after using bridge loan to cover debt, actually end
% up saving in net.

ar_b_bridge = zeros(size(ar_coh_today));

if (bl_b_is_principle)

    % c_bridge is cost of borrowing in next period consumption
    ar_b_bridge(ar_coh_today<0) = ar_coh_today(ar_coh_today<0);
    ar_c_bridge = ar_b_bridge.*(1+fl_r_bridge);

else

    % c_bridge is the gain from borrowing in current period consumption
    ar_b_bridge(ar_coh_today<0) = ar_coh_today(ar_coh_today<0).*(1+fl_r_bridge);
    ar_c_bridge = (-1)*ar_b_bridge./(1+fl_r_bridge);

end

% remaining aprime after allocating to pay debt not covered by coh
ar_aprime_nobridge = ar_aprime - ar_b_bridge;

if (bl_display_infbridge)
    tab_aprime_bridge = table(ar_coh_today, ar_b_bridge, ar_c_bridge, ar_aprime, ar_aprime_nobridge);
    disp('coh_today: includes all income - debt');
    disp(tab_aprime_bridge);
end

% Display
if (bl_display_infbridge)
    tab_bridge = table(ar_coh_today, ar_aprime, ar_aprime_nobridge, ar_c_bridge, ar_b_bridge);
    disp(tab_bridge);
end

end
