%% 
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository> 
% Table of Content.*

function [ar_a] = ffs_az_gen_borrsave_grid(varargin)
%% FFS_AZ_GEN_BORRSAVE_GRID get funcs, params, states choices shocks grids
% Generate savings and borrowing states/choice grids. There are three types
% of grids:
% # savings choice/state grid
% # borrowing choice/state grid when default is not allowed
% # borrowing choice/state grid when default is allowed
%
% The parameters for this function and the role of cmin, ymin, default, etc
% are discussed on this page:
% <https://fanwangecon.github.io/CodeDynaAsset/docs/README_cminymin_borrsave.html
% cmin ymin borr save>.
%
% @param fl_b_bd float exogenously set borrowing bound. If zero, borrowing
% is not allowed. If lower than zero, means borrowing is allowed. When
% borrowing is allowed fl_b_bd is only the actual bound on borrowing if it
% is tighter than the natural borrowing constraint based on minimum income
% in the case where defaults are not allowed. And if it is tighther than
% the defaults borrowing constraint based on maximium income in the case
% where defaults are allowed.  
%
% @param bl_default boolean if fl_b_bd is below zero, then is defaults
% allowed or not? What does default mean? No-default means households under
% all shocks, can always repay existing debts and interests with new debts.
% This means banks do not have to worry about non-repayment. Households
% always pay back what they are owed. Note that debt roll-over is allowed
% here. What does default mean, it means under some states of shocks,
% households face debt that they can not repay, so they go to minimum
% consumption, which is utility for default. Their optimal choice is
% savings = 0 for the next period. 
%
% @param ar_z array array of exogenous income shocks for the exogeous
% shocks to inelastic labor supply 
%
% @return fl_w float wage
%
% @return fl_a_max float maximum savings level
%
% @return it_a_n integer number of save/borrow grid points
%
% @example
%
% @include
%

cl_params_len = length(varargin);

%% Default Folder Parameters

fl_b_bd = -100;
bl_default = 1;
ar_z = [0.3, 2.3];
fl_w = 1;
fl_r_borr = 0.05;
fl_a_max = 50;
it_a_n = 100;

cl_params = {fl_b_bd bl_default ar_z ...
                   fl_w fl_r_borr fl_a_max it_a_n};

%% Parse Parameters
% numvarargs is the number of varagin inputted
[cl_params{1:cl_params_len}] = varargin{:};
fl_b_bd = cl_params{1};
bl_default = cl_params{2};
ar_z = cl_params{3};
fl_w = cl_params{4};
fl_r_borr = cl_params{5};
fl_a_max = cl_params{6};
it_a_n = cl_params{7};

%% Savings Only
if (fl_b_bd >= 0)
    fl_a_min = 0;
end

%% Borrowing
if (fl_b_bd < 0)
    
    %% Borrowing not allowing for default
    if (~bl_default)
        fl_ar_z_min = min(ar_z);
        fl_y_min = fl_ar_z_min*fl_w;
        fl_a_min = max(fl_b_bd, -fl_y_min/fl_r_borr);
    end
    
    %% Borrowing allowing for default
    if (bl_default)
        fl_ar_z_max = max(ar_z);
        fl_y_max = fl_ar_z_max*fl_w;
        fl_a_min = max(fl_b_bd, -fl_y_max/fl_r_borr);
    end
    
    
end

%% Grid
ar_a = linspace(fl_a_min, fl_a_max, it_a_n);

%% Add Zero
ar_a = [0 ar_a];
ar_a = sort(unique(ar_a));

end
