%%
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

function [mt_coh] = ffs_fibs_coh(varargin)
%% FFS_FIBS_CO
% Suppose I want to borrow 5 dollar today. That will increase current
% consumption by 5 dollars, but its cost to next period cash-on-hand is a
% function of how I get the 5 dollar of borrowing, whether I borrow from
% formal sources, informal or joint sources. Suppose the formal borrowing
% grid has three points, you can borrow 0, 3 or 6, these are the options:
%
% # borrow all 5 informally
% # borrow 3 formally and 2 informally
% # borrow 6 formally and save 1
% # additionally, if 5 is actually on the grid, you would borrow 5 only
% formally if that offers the lower interest rate. 
% 
% This function here finds which one of of these options is optimal
% conditional on aprime_nobridge choice. Households choose among the different
% aprime_nobridge choices given these within aprime_nobridge choice optimal formal/informal
% allocations. 
%
% *aprime_nobridge* is the borrowing choices that do not include bridge
% loan. The additional borrowing on top of bridge loans.
%
% @param fl_r_inf float informal interest rate
%
% @param fl_r_fsv float (formal) savings interest rate
%
% @param fl_r_fbr float borrowing interest rate
%
% @param aprime_nobridge array level of aggregate borrowing excluding
% bridge loan. Note that bridge loan is needed if coh is negative and
% households can not pay back principle and interests.
%
% @return max_c_nextperiod array maximum next period consumption given aprime_nobridge and
% optimal formal informal choices
%
% @return fl_for_borr array formal borrowing choices that maximized next
% period consumption for each element of aprime_nobridge vector
%
% @return fl_for_save array (formal) savings choices that maximized next
% period consumption for each element of aprime_nobridge vector
%
% @return fl_inf_borr array informal borrowing choices that maximized next
% period consumption for each element of aprime_nobridge vector
%
% @example
%
%   [max_c, fl_for_borr, fl_for_save, fl_inf_borr, fl_coh_add] = ...
%        ffs_akz_set_functions(fl_r_inf, fl_r_fsv, fl_r_fbr, fl_for_br_block, aprime_nobridge);
%

end