%% Generating model functions
function [f_util_log, f_util_crra, f_util_standin, f_inc, f_coh, f_cons] = ffs_az_set_functions(varargin)
% @example
% [f_util, fl_u_neg_c, f_prod, f_coh, f_cons, f_v] = equa_owkbz(fl_crra, fl_u_neg_c, ...
%            fl_Amean, fl_alpha, fl_delta, fl_r, ...
%            fl_v_alpha, fl_v_beta, fl_b_bd)

%% Catch Error
optional_params_len = length(varargin);
if optional_params_len > 10
    error('grid_wkb:TooManyOptionalParameters', ...
        'allows at most 10 optional parameters');
end

%% Default Parameters
[fl_crra, fl_c_min] = deal(1, 0.001);
[fl_r_save, fl_r_borr, fl_wage] = deal(0.02, 0.02, 1.23);
optional_params = {fl_crra fl_c_min fl_r_save fl_r_borr fl_wage};

%% Parse Parameters
% numvarargs is the number of varagin inputted
[optional_params{1:optional_params_len}] = varargin{:};
[fl_crra, fl_c_min, fl_r_save, fl_r_borr, fl_wage] = optional_params{:};

%% Equations
% utility
f_util_log = @(c) log(c);
f_util_crra = @(c) (((c).^(1-fl_crra)-1)./(1-fl_crra));
% Production Function
f_inc = @(z, b) (z*fl_wage + b.*(fl_r_save).*(b>0) + b.*(fl_r_borr).*(b<=0)); % z already exp
% Cash on Hand, b is principle and interest
f_coh = @(z, b) (z*fl_wage + b.*(1+fl_r_save).*(b>0) + b.*(1+fl_r_borr).*(b<=0));
% Simple Consumption b and k
f_cons = @(z, b, bprime) (f_coh(z, b) - bprime);
% Simple Consumption b and k
f_util_standin = @(z, b) f_util_log(f_coh(z,b).*(f_coh(z,b) > 0) + ...
                                    fl_c_min.*(f_coh(z,b) <= 0));

end
