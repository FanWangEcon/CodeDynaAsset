%% Set Model Functions
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [f_util_log, f_util_crra, f_util_standin, f_inc, f_coh, f_cons] = ffs_az_set_functions(varargin)
%% FFS_AZ_SET_FUNCTIONS setting model functions
% define functions here to avoid copy paste mistakes
%
% @param fl_crra float crra utility
%
% @param fl_r_save float savings interest rate
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
%        ffs_az_set_functions(fl_crra, fl_r_save, fl_w);
%
% @seealso
%
% * initialize parameters: <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html ffs_az_set_default_param>
% * initialize functions: <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_functions.html ffs_az_set_functions>
% * set shock wage: <https://github.com/FanWangEcon/CodeDynaAsset/blob/master/tools/ffto_gen_tauchen_jhl.m ffto_gen_tauchen_jhl>
% * gateway function processing grid, paramters, functions: <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_get_funcgrid.html ffs_az_get_funcgrid>
%

%% Default

[fl_crra] = deal(1);
[fl_r_save, fl_w] = deal(0.02, 1.23);
default_params = {fl_crra fl_r_save fl_w};

%% Parse Parameters

% numvarargs is the number of varagin inputted
[default_params{1:length(varargin)}] = varargin{:};
[fl_crra, fl_r_save, fl_w] = default_params{:};

%% Equations

% utility
f_util_log = @(c) log(c);
f_util_crra = @(c) (((c).^(1-fl_crra)-1)./(1-fl_crra));
% Production Function
f_inc = @(z, b) (z*fl_w + b.*(fl_r_save)); % z already exp
% Cash on Hand, b is principle
f_coh = @(z, b) (z*fl_w + b.*(1+fl_r_save));
% Simple Consumption b and k
f_cons = @(z, b, bprime) (f_coh(z, b) - bprime);
% Simple Consumption b and k
f_util_standin = @(z, b) f_util_log(f_coh(z,b));

end
