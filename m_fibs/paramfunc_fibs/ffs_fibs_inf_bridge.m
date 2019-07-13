%% Formal Informal Bridge Loan Handling
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
%%

function [ar_aprime_nobridge, ar_b_bridge, ar_c_bridge] = ffs_fibs_inf_bridge(varargin)
%% FFS_FIBS_INF_BRIDGE Amount of Informal Borrowing Needed as Bridge Loans
% Bridge loan needed to pay for debt that is still unpaid due to
% insufficient cash on hand. Potentially, only informal lender is willing
% to extend this loan offers. So some existing debts are paid back with
% revenue, parts that revenue can not cover is paid back potentially with
% informal borrowing. This works with single and multiple assets.
%
% @param bl_b_is_principle boolean solving with aggregate save/borr as
% principles + interests, or just principles no interests. If true,
% principels only, no interests. These refer to the _ar_aprime_ vector.
%
% @param fl_r_bridge float interest rate for bridge loan
%
% @param fl_r_fsv float (formal) savings interest rate
%
% @param fl_r_fbr float borrowing interest rate
%
% @param ar_aprime array N by 1 level of aggregate borrowing including
% potentially bridge loan. Note that bridge loan is needed if coh is
% negative and households can not pay back principle and interests.
%
% @param ar_coh_today array N by 1 the level of cash-on-hand today, when the
% borrowing and savings decisions are made. If this is positive, then
% households freely borrow, do not need bridge loans. If this is negative
% households need to first borrow to meet bridge loan needs.
%
% @return ar_aprime_nobridge array next period asset choice without debt
% incurred for bridge loans. This is the difference between _ar_coh_today_
% and _ar_b_bridge_.
%
% @return ar_b_bridge array bridge loan debt to pay for unpaid uncovered
% cash-on-hand. Includes interests and principle if _bl_b_is_principle_ is
% false. If _bl_b_is_principle_ is true, only includes principle.
%
% @return ar_c_bridge array consumption gain (for today) from bridge loan,
% or consumption costs (for next period) of bridge loans. Both consider the
% impact of principles as well as interest rates on bridge loans.
%
% @example
%
%   bl_input_override = true;
%   [ar_aprime_nobridge, ar_c_bridge] = ...
%        ffs_fibs_inf_bridge(bl_b_is_principle, fl_r_bridge, ...
%                              ar_aprime, ar_coh_today, ...
%                              bl_display_infbridge, bl_input_override);
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

if (~isempty(varargin))
    
    % override when called from outside
    [bl_b_is_principle, fl_r_bridge, ar_aprime, ar_coh_today, bl_display_infbridge] = varargin{:};
    
else
    close all

    % Default
    it_param_set = 4;
    [param_map, ~] = ffs_ipwkbz_fibs_set_default_param(it_param_set);

    % Gather Inputs from param_map
    params_group = values(param_map, {'bl_b_is_principle', 'fl_r_inf'});
    [bl_b_is_principle, fl_r_inf] = params_group{:};
    % bl_b_is_principle = true;
    
    % For benchmark, assume that the informal lender
    fl_r_bridge = fl_r_inf;

    % Testing COH and Aprime Vectors
    ar_aprime =    [-20. -5,-5, -4.5,-4.5, -0.1,-0.1]';
    ar_coh_today = [-19, 1, -1,   1,-1,      1,0 ]';
    
%     % Testing COH and Aprime coh is row vector, aprime is column vector
%     ar_aprime =    linspace(-5, 10, 5)';
%     ar_coh_today = linspace(-10, 0, 10);

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

%% Display
if (bl_display_infbridge)
        
    disp(['bl_b_is_principle:', num2str(bl_b_is_principle)]);
    disp(['fl_r_bridge:', num2str(fl_r_bridge)]);
    disp(['ar_aprime:', num2str(ar_aprime')]);
    disp(['ar_b_bridge:', num2str(ar_coh_today')]);
    disp(['bl_display_infbridge:', num2str(bl_display_infbridge)]);

    tab_aprime_bridge = table(ar_coh_today, ar_aprime, ar_b_bridge, ar_c_bridge, ar_aprime_nobridge);
    tab_aprime_bridge.Properties.VariableDescriptions{'ar_coh_today'} = ...
        '*ar_coh_today*: cash on hand someone arrives in the period with given debt and current income';
    tab_aprime_bridge.Properties.VariableDescriptions{'ar_aprime'} = ...
        '*ar_aprime*: func called during finding optimal aprime, this is the current aprime overall choice';
    tab_aprime_bridge.Properties.VariableDescriptions{'ar_b_bridge'} = ...
        '*ar_b_bridge*: amount of bridge loan required to cover negative coh (includes interest if bl_b_is_principle = false)';
    tab_aprime_bridge.Properties.VariableDescriptions{'ar_c_bridge'} = ...
        '*ar_c_bridge*: consumption gain today from the bridge loan to cover negative coh; or consumption cost tomorrow for debt which increases c today';
    tab_aprime_bridge.Properties.VariableDescriptions{'ar_aprime_nobridge'} = ...
        ['*ar_aprime_nobridge*ar_aprime_nobridge:' ...
         'aprime = -10, -5 for bridge; -5 left for other borrowing choices;' ...
         'aprime = -10, -11 for bridge (given r), +1 savings left, reduces consumption, back to neg coh, infeasible state;'];
    
    disp(tab_aprime_bridge);
    sc_summary = summary(tab_aprime_bridge);
    cl_var_name = fieldnames(sc_summary);
    for it_var_name = 1:length(cl_var_name)
        disp(sc_summary.(cl_var_name{it_var_name}).Description);
    end
    
end

end
