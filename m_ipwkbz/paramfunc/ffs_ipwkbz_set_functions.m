%% Set Model Functions (Interpolated + Percentage + Risky + Safe Asset + Save + Borrow)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [f_util_log, f_util_crra, f_util_standin, f_prod, f_inc, f_coh, f_cons] = ffs_ipwkbz_set_functions(varargin)
%% FFS_IPWKBZ_SET_FUNCTIONS setting model functions
% define functions here to avoid copy paste mistakes. This function is
% identical to the
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_functions.html
% ffs_akz_set_functions> file.
%
% @param fl_crra float crra utility
%
% @param fl_c_min float minimum consumption
%
% @param fl_Amean float mean productivity level
%
% @param fl_alpha float Cobb-Douglas elasticity crs if = 1 and is risky
% financial investment like stocks, if decreasing return to scales, could
% be interpreted as risky physical capital investment.
%
% @param fl_delta float risky capital depreciation if this is risky
% financial investment, deprecitation could be full, there is no additional
% return from risky financial investment except for its production function
% earnings. If we are looking are risky capital investment, there is a
% production function return, output. Then we can also sell the risky
% capital with some depreciation. fl_delta = 0.1 means 10 percent
% depreciation.
%
% @param fl_r_save float savings interest rate
%
% @param fl_r_borr float borrowing interest rate
%
% @param fl_w float wage rate
%
% @return f_util_log handle log utility
%
% @return f_util_crra handle crra utility general
%
% @return f_util_standin handle log utility with coh for testing graphing codes
%
% @return f_inc income handle equation wage and interests
%
% @return f_coh handle cash on hand equation given current period shock
%
% @return f_cons handle consumption equation given coh and asset choice
%
% @example
%
%   [f_util_log, f_util_crra, f_util_standin, f_inc, f_coh, f_cons] = ...
%        ffs_akz_set_functions(fl_crra, fl_c_min, fl_Amean, fl_alpha, fl_delta, ...
%                              fl_r_save, fl_r_borr, fl_w);
%

%% Default

[fl_crra, fl_c_min, fl_b_bd] = deal(1.5, 0.001, -20);
[fl_Amean, fl_alpha, fl_delta] = deal(1, 0.36, 0.08);
[fl_r_save, fl_r_borr, fl_w] = deal(0.02, 0.02, 1.23);
default_params = {fl_crra fl_c_min fl_b_bd ...
    fl_Amean fl_alpha fl_delta fl_r_save fl_r_borr fl_w};


%% Parse Parameters

% numvarargs is the number of varagin inputted
[default_params{1:length(varargin)}] = varargin{:};
[fl_crra, fl_c_min, fl_b_bd,...
    fl_Amean, fl_alpha, fl_delta, fl_r_save, fl_r_borr, fl_w] = default_params{:};

%% Equations Utility

f_util_log = @(c) log(c);
f_util_crra = @(c) (((c).^(1-fl_crra)-1)./(1-fl_crra));

%% Equations Production
% production function, z already exp, possible decreasing return to scale.
% if fl_alpha = 1, crs. That means we have a risky asset like stocks vs
% safe asset like bond. If we have decrease return to scale, can be
% interpreted as capital investments.

f_prod = @(z, k) ((fl_Amean.*(z)).*(k.^(fl_alpha)));

%% Equations Income
% Income equation now based on b = interest + principles, fl_w is fixed
% wage, shock now on risky capital. Three sources of income, production
% income minus depreciation, wage income (could be zero), also interest
% income or interest costs dependong on borrowing or savings

% f_inc = @(z, b, k) (f_prod(z, k) - (fl_delta)*k ...
%                     + fl_w ...
%                     + (b./(1+fl_r_save)).*fl_r_save.*(b>0) + (b/(1+fl_r_borr)).*(fl_r_borr).*(b<=0)); % z already exp

f_inc = @(z, b, k) (f_prod(z, k) - (fl_delta)*k ...
                     + fl_w ...
                     + b.*(fl_r_save).*(b>0) + b.*(fl_r_borr).*(b<=0)); % z already exp

%% Equations Cash-on-Hand
% Cash on Hand, b is principle and interest. Very important to include fl_w in
% the COH equation, with fl_w > 0, there will be an income floor, lowest coh
% that we will solve for will not be 0. with fl_w > 0, the c_min parameter is
% only there to reset invalid choice grid values, have no effects on value func.
% f_coh = @(z, b, k) (f_prod(z, k) + k*(1-fl_delta) + fl_w + b);
f_coh = @(z, b, k) (f_prod(z, k) + k*(1-fl_delta) + fl_w + b.*(1+fl_r_save).*(b>0) + b.*(1+fl_r_borr).*(b<=0));

%% Equations Consumption
% Simple Consumption given cash-on-hand
% f_cons = @(coh, bprime, kprime) (coh - kprime - ((bprime./(1+fl_r_save)).*(bprime>0)) - ((bprime./(1+fl_r_borr)).*(bprime<=0)));
f_cons = @(coh, bprime, kprime) (coh - kprime - bprime);

%% Equations Stand-in Fake Utility for Graphs
% Utility for graphing with random data, note that when we graph with coh
% as the state variable using this equation here, there is no effect of
% shock on utility, it is fully captured by the coh.
f_util_standin = @(z, b, k) f_util_log((f_coh(z,b,k)-fl_b_bd).*((f_coh(z,b,k) - fl_b_bd) > fl_c_min) + ...
                                       fl_c_min.*((f_coh(z,b,k) - fl_b_bd) <= fl_c_min));


end
