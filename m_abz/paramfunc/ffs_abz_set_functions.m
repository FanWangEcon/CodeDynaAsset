%% Set Model Functions One Asset Borrowing and Savings
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [f_util_log, f_util_crra, f_util_standin, f_awithr_to_anor, f_coh, f_cons_coh, f_cons_checkcmin] = ...
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
%   [f_util_log, f_util_crra, f_util_standin, f_inc, f_coh, f_cons_coh, f_cons, f_cons_checkcmin] = ...
%        ffs_abz_set_functions(fl_crra, fl_c_min, fl_r_save, fl_w);
%
% @seealso
%
% * initialize parameters: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_set_default_param.html ffs_abz_set_default_param>
% * initialize functions: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_set_functions.html ffs_abz_set_functions>
% * set asset grid: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_gen_borrsave_grid.html ffs_abz_gen_borrsave_grid>
% * set shock borrow rate: <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_gen_discrete_var.html fft_gen_discrete_var>
% * set shock wage: <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/tools/ffto_gen_tauchen_jhl.m ffto_gen_tauchen_jhl>
% * gateway function processing grid, paramters, functions: <https://fanwangecon.github.io/CodeDynaAsset/m_abz/paramfunc/html/ffs_abz_get_funcgrid.html ffs_abz_get_funcgrid>
%

%% Default

[fl_crra, fl_c_min] = deal(1, 0.001);
[fl_r_save, fl_w] = deal(0.02, 1.23);
default_params = {fl_crra fl_c_min fl_r_save fl_w};

%% Parse Parameters

% numvarargs is the number of varagin inputted
[default_params{1:length(varargin)}] = varargin{:};
[fl_crra, fl_c_min, fl_r_save, fl_w] = default_params{:};

%% Equations if Choice is Principle
% Choice is principle does not work in this structure with variation in
% borrowing interest rate. 

% utility
f_util_log = @(c) log(c);
f_util_crra = @(c) (((c).^(1-fl_crra)-1)./(1-fl_crra));

% COH and C
f_coh = @(z, b) (z*fl_w + b);
f_awithr_to_anor = @(fl_r_borr, bprime) (bprime.*(1./(1+fl_r_save)).*(bprime>0) + bprime.*(1./(1+fl_r_borr)).*(bprime<=0));
f_cons_coh = @(coh, fl_r_borr, bprime) (coh - f_awithr_to_anor(fl_r_borr, bprime));


% Support
f_coh_princ = @(fl_r_borr, z, b) (z*fl_w + (b.*(1+fl_r_save).*(b>0) + b.*(1+fl_r_borr).*(b<=0)));
f_cons_coh_princ = @(coh, bprime) (coh - bprime);

% f_cons_checkcmin is not used in main solution, but at the end of solution file to
% get consumption based on optimal choices. Here we add conditioning based
% on cmin, so that consumption is not negative. If do not do this, for
% defaulters, next period a'=0, would seem like they have large negative
% consumption.
f_cons_checkcmin = @(fl_r_borr, z, b, bprime) ((f_cons_coh(f_coh(z, b),fl_r_borr,bprime)).*((f_cons_coh(f_coh(z, b),fl_r_borr,bprime)) >= fl_c_min) + ...
                                               fl_c_min.*((f_cons_coh(f_coh(z, b),fl_r_borr,bprime)) < fl_c_min));
% Simple Consumption b and k
f_util_standin = @(fl_r_borr, z, b) f_util_log(f_coh_princ(fl_r_borr,z,b).*(f_coh_princ(fl_r_borr,z,b) > 0) + ...
                                    fl_c_min.*(f_coh_princ(fl_r_borr,z,b) <= 0));


end
