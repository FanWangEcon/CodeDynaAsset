%% Set Model Functions
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [f_util_log, f_util_crra, f_util_standin, f_inc, f_coh, f_cons_coh, f_cons, f_cons_checkcmin] = ...
    ffs_abz_set_functions(varargin)
%% FFS_ABZ_SET_FUNCTIONS setting model functions
% define functions here to avoid copy paste mistakes
%
% @param fl_crra float crra utility
%
% @param fl_c_min float minimum consumption
%
% @param fl_r_save float savings interest rate, the borrowing interest rate
% is a paramter for the coh and inc functions, not a parameter for the
% overall function here.
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
%        ffs_abz_set_functions(fl_crra, fl_c_min, fl_r_save, fl_r_borr, fl_w);
%

%% Default

[fl_crra, fl_c_min] = deal(1, 0.001);
[fl_r_save, fl_w] = deal(0.02, 1.23);
default_params = {fl_crra fl_c_min fl_r_save fl_w};

%% Parse Parameters

% numvarargs is the number of varagin inputted
[default_params{1:length(varargin)}] = varargin{:};
[fl_crra, fl_c_min, fl_r_save, fl_w] = default_params{:};

%% Equations
% z here refers to z_wage

% utility
f_util_log = @(c) log(c);
f_util_crra = @(c) (((c).^(1-fl_crra)-1)./(1-fl_crra));
% Production Function
f_inc = @(fl_r_borr, z, b) (z*fl_w + (b.*(fl_r_save).*(b>0) + b.*(fl_r_borr).*(b<=0))); % z already exp
% Cash on Hand, b is principle
f_coh = @(fl_r_borr, z, b) (z*fl_w + (b.*(1+fl_r_save).*(b>0) + b.*(1+fl_r_borr).*(b<=0)));
% Simple Consumption b and k
f_cons_coh = @(coh, bprime) (coh - bprime);
f_cons = @(z, b, bprime) (f_coh(z, b) - bprime);
% f_cons_checkcmin is not used in main solution, but at the end of solution file to
% get consumption based on optimal choices. Here we add conditioning based
% on cmin, so that consumption is not negative. If do not do this, for
% defaulters, next period a'=0, would seem like they have large negative
% consumption.
f_cons_checkcmin = @(fl_r_borr, z, b, bprime) ((f_coh(fl_r_borr, z, b) - bprime).*(f_coh(fl_r_borr, z, b) - bprime >= fl_c_min) + ...
                                               fl_c_min.*(f_coh(fl_r_borr, z, b) - bprime < fl_c_min));
% Simple Consumption b and k
f_util_standin = @(fl_r_borr, z, b) f_util_log(f_coh(fl_r_borr,z,b).*(f_coh(fl_r_borr,z,b) > 0) + ...
                                    fl_c_min.*(f_coh(fl_r_borr,z,b) <= 0));

end
