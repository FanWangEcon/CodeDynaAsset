%% Set Model Functions (Interpolated + Percentage + Risky + Safe Asset + FIBS + RShock)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [f_util_log, f_util_crra, f_util_standin, f_util_standin_coh, ...
    f_prod, f_inc, f_coh, f_coh_fbis, f_coh_save, f_cons] = ffs_ipwkbzr_fibs_set_functions(varargin)
%% FFS_IPWKZR_FIBS_SET_FUNCTIONS setting model functions

%% Default

[fl_crra, fl_c_min, fl_b_bd] = deal(1.5, 0.001, -20);
[fl_Amean, fl_alpha, fl_delta] = deal(1, 0.36, 0.08);
[fl_w, fl_r_fbr, fl_r_fsv] = deal(1.28, 0.035, 0.025);
bl_display_setfunc = false;

default_params = {fl_crra fl_c_min fl_b_bd ...
    fl_Amean fl_alpha fl_delta fl_w fl_r_fbr fl_r_fsv bl_display_setfunc};

%% Parse Parameters

% numvarargs is the number of varagin inputted
[default_params{1:length(varargin)}] = varargin{:};
[fl_crra, fl_c_min, fl_b_bd, fl_Amean, fl_alpha, fl_delta, ...
    fl_w, fl_r_fbr, fl_r_fsv, bl_display_setfunc] = default_params{:};

%% Equations Utility

f_util_log = @(c) log(c);
f_util_crra = @(c) (((c).^(1-fl_crra)-1)./(1-fl_crra));

%% Equations Production

f_prod = @(z, k) ((fl_Amean.*(z)).*(k.^(fl_alpha)));

%% Equations Income

f_inc = @(z, k, fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save) ...
    (f_prod(z, k) - (fl_delta)*k + fl_w + ...
        (ar_for_borr.*(fl_r_fbr) + ar_inf_borr.*(fl_r_inf) + ...
         ar_for_save.*(fl_r_fsv)));

%% Equations Cash-on-Hand
%
% *b_cons_cost*: For borrowing, given the overall b choice, function
% <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_min_c_cost_bridge.html
% ffs_fibs_min_c_cost_bridge> finds the optimal formal and informal joint
% choices, and we obtain the consumption cost of borrowing from that
%
% Unlike for the *abz (fibs)* problem, where the equations below are in
% terms of consumption, here we know b is principle due to two stage
% solution method.
%

% coh, where the b_with_r include various formal and informal choices with
% interest rate.
f_coh = @(z, b_with_r, k) (f_prod(z, k) + k*(1-fl_delta) + fl_w + b_with_r);

% If borrowing overall, various files in
f_coh_fbis = @(fl_r_inf, ar_for_borr, ar_inf_borr, ar_for_save) ...
                (ar_for_borr.*(1+fl_r_fbr) + ...
                 ar_inf_borr.*(1+fl_r_inf) + ...
                 ar_for_save.*(1+fl_r_fsv));

% If saving overall
f_coh_save = @(b) (b.*(1+fl_r_fsv));

%% Equations Consumption (b=principle)
% Note this is different from
% <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/m_abz_paramfunc/html/ffs_abz_fibs_set_functions.html
% ffs_abz_fibs_set_functions>, here only need to deal with b = principle,
% otherwise two stage solution does not work. So from today's perspective,
% choosing k and b for next period in principles.

f_cons = @(coh, bprime, kprime) (coh - kprime - bprime);

%% Equations Stand-in Fake Utility for Graphs
% Utility for graphing with random data, note that when we graph with coh
% as the state variable using this equation here, there is no effect of
% shock on utility, it is fully captured by the coh.

f_util_standin = @(z, b, k) f_util_log((f_coh(z,b,k)-fl_b_bd).*((f_coh(z,b,k) - fl_b_bd) > fl_c_min) + ...
                                        fl_c_min.*((f_coh(z,b,k) - fl_b_bd) <= fl_c_min));

f_util_standin_coh = @(coh, fl_r_borr) f_util_log((coh-fl_b_bd).*( (coh > 0) & (((coh - fl_b_bd)./(1)) > fl_c_min)) + ...
                                                  ((coh-fl_b_bd)./(1)).*( (coh <= 0) & (((coh - fl_b_bd)./(1)) > fl_c_min)) + ...
                                                  (fl_c_min./(1+fl_r_borr)).*( (((coh - fl_b_bd)./(1)) <= fl_c_min)));


end
