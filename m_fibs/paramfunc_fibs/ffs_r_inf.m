%%
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

function [mt_r_inf_trans, ar_r_inf_stationary, ar_r_inf] = ffs_r_inf(varargin)
%% FFS_R_INF informal interest rate distribution
% Households draw informal interest rates. Probability for each discrete
% element in the sample space of the discrete random variable for informal
% interest rate could be estimated.
%
% @param fl_r_inf float informal interest rate
%
% @param fl_r_fsv float (formal) savings interest rate
%
% @param fl_r_fbr float borrowing interest rate
%
% @param aprime array level of aggregate borrowing excluding
% bridge loan. Note that bridge loan is needed if coh is negative and
% households can not pay back principle and interests.
%
% @param coh_today array the level of cash-on-hand today, when the
% borrowing and savings decisions are made. If this is positive, then
% households freely borrow, do not need bridge loans. If this is negative
% households need to first borrow to meet bridge loan needs. 
%
% @return mt_r_inf_trans transition matrix for the informal interest rate
% discrete random variable
%
% @return ar_r_inf_stationary stationary distribution of the informal
% interest rate discrete random variable
%
% @return ar_r_inf values associated with each element in the sample space
% of the informal interest rate discrete random variable
%
% @example
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
