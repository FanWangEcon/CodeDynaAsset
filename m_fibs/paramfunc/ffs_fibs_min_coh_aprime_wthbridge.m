%%
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

function [max_c_nextperiod_nobridge, fl_for_borr, fl_for_save, fl_inf_borr] = ffs_fibs_min_coh_aprime_wthbridge(varargin)
%% FFS_FIB_MIN_COH_aprime optimal borrow for+inf given borrowing
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
% conditional on aprime choice. Households choose among the different
% aprime choices given these within aprime choice optimal formal/informal
% allocations. 
%
% *aprime* is the borrowing choices that do not include bridge
% loan. The additional borrowing on top of bridge loans.
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
% @return max_c_nextperiod array maximum next period consumption given aprime and
% optimal formal informal choices
%
% @return fl_for_borr array formal borrowing choices that maximized next
% period consumption for each element of aprime vector
%
% @return fl_for_save array (formal) savings choices that maximized next
% period consumption for each element of aprime vector
%
% @return fl_inf_borr array informal borrowing choices that maximized next
% period consumption for each element of aprime vector
%
% @example
%
%   [max_c, fl_for_borr, fl_for_save, fl_inf_borr, fl_coh_add] = ...
%        ffs_akz_set_functions(fl_r_inf, fl_r_fsv, fl_r_fbr, fl_for_br_block, aprime);
%

%% Default
[fl_r_inf, fl_r_fsv, fl_r_fbr] = deal(0.15, 0.01, 0.06);
[fl_for_br_block] = -1;
[aprime] =    [-5,-5, -4.5,-4.5, -0.1,-0.1];
[coh_today] = [ 1, -1,   1,-1,      1,0 ];
bl_display_fibs_min_coh_aprime = true;
default_params = {fl_r_inf fl_r_fsv fl_r_fbr fl_for_br_block aprime coh_today bl_display_fibs_min_coh_aprime};

%% Parse Parameters

% numvarargs is the number of varagin inputted
[default_params{1:length(varargin)}] = varargin{:};
[fl_r_inf, fl_r_fsv, fl_r_fbr, fl_for_br_block, aprime, coh_today, ...
    bl_display_fibs_min_coh_aprime] = default_params{:};

%% Bridge Loan Required

b_inf_bridge = zeros(size(coh_today));
b_inf_bridge(coh_today<0) = coh_today(coh_today<0);
c_bridge = b_inf_bridge.*(1+fl_r_inf);
aprime_nobridge = aprime - b_inf_bridge;

if (bl_display_fibs_min_coh_aprime)
    tab_aprime_bridge = table(coh_today', b_inf_bridge', c_bridge', aprime', aprime_nobridge');
    tab_aprime_bridge.Properties.VariableNames = {'coh_today', 'b_inf_bridge', 'c_bridge', 'aprime', 'aprime_nobridge'};
    disp('coh_today: includes all income - debt');
    disp(tab_aprime_bridge);    
end

%% Compute Consumption Informal Borrowing only

b_infonly_inf = (aprime_nobridge);
c_infonly = b_infonly_inf.*(1+fl_r_inf);

if (bl_display_fibs_min_coh_aprime)
    tab_c_infonly = table(aprime_nobridge', b_infonly_inf', c_infonly');
    tab_c_infonly.Properties.VariableNames = {'aprime_nobridge','b_infonly_inf','c_infonly'};
    disp(['informal borrow interest (fl_r_inf):', num2str(fl_r_inf)]);
    disp(tab_c_infonly);    
end

%% Compute Consumption Values Formal + Informal Borrowing Jointly
% if right on grid, informal could be zero. 

for_br_block_dvd = aprime_nobridge./fl_for_br_block;
for_br_block_dvd_floor = abs(floor(for_br_block_dvd).*fl_for_br_block);
for_br_block_dvd_ceil = abs(ceil(for_br_block_dvd).*fl_for_br_block);

b_infforb_inf = (aprime_nobridge + for_br_block_dvd_floor);
b_infforb_for = -(for_br_block_dvd_floor);
c_infforb =  (b_infforb_inf.*(1+fl_r_inf) + b_infforb_for.*(1+fl_r_fbr));

if (bl_display_fibs_min_coh_aprime)
    tab_c_infforb = table(aprime_nobridge', b_infforb_for', b_infforb_inf', c_infforb');
    tab_c_infforb.Properties.VariableNames = {'aprime_nobridge','b_infforb_for','b_infforb_inf','c_infforb'};
    disp(['formal borrow interest (fl_r_fbr):', num2str(fl_r_fbr)]);
    disp(['informal borrow interest (fl_r_inf):', num2str(fl_r_inf)]);
    disp(['formal block size (fl_for_br_block):', num2str(fl_for_br_block)]);
    disp(tab_c_infforb);    
end

%% Compute Consumption Values Formal Borrowing + Formal Savings

b_forbrsv_sav = (aprime_nobridge + for_br_block_dvd_ceil);
b_forbrsv_brr = -(for_br_block_dvd_ceil);
c_forbrsv = (b_forbrsv_sav*(1+fl_r_fsv) + b_forbrsv_brr*(1+fl_r_fbr));

if (bl_display_fibs_min_coh_aprime)
    tab_c_forbrsv = table(aprime_nobridge', b_forbrsv_brr', b_forbrsv_sav', c_forbrsv');
    tab_c_forbrsv.Properties.VariableNames = {'aprime_nobridge','b_forbrsv_brr','b_forbrsv_sav','c_infforb'};
    disp(['formal borrow interest (fl_r_fbr):', num2str(fl_r_fbr)]);
    disp(['savings interest (fl_r_fsv):', num2str(fl_r_fsv)]);
    disp(['formal block size (fl_for_br_block):', num2str(fl_for_br_block)]);
    disp(tab_c_forbrsv);    
end

%% Maximize Consumption For non-Bridge Loan Component
% Trying to minimize consumption costs
[max_c_nextperiod_nobridge, max_idx] = max([c_infonly; c_infforb; c_forbrsv],[], 1);

if (bl_display_fibs_min_coh_aprime)
    tab_c_max = table(aprime_nobridge', max_c_nextperiod_nobridge', max_idx', c_infonly', c_infforb', c_forbrsv');
    tab_c_max.Properties.VariableNames = {'aprime_nobridge', 'max_c_nextperiod_nobridge', 'max_idx', 'c_infonly','c_infforb','c_forbrsv'};
    disp(['informal borrow interest (fl_r_inf):', num2str(fl_r_inf)]);    
    disp(['formal borrow interest (fl_r_fbr):', num2str(fl_r_fbr)]);
    disp(['savings interest (fl_r_fsv):', num2str(fl_r_fsv)]);
    disp(['formal block size (fl_for_br_block):', num2str(fl_for_br_block)]);
    disp(tab_c_max);    
end

%% Consumption Cost of Aprime Choice Bridge Loan includes
aprime_c_cost_total = c_bridge + max_c_nextperiod_nobridge;

if (bl_display_fibs_min_coh_aprime)
    tab_c_max = table(coh_today', aprime', aprime_nobridge', aprime_c_cost_total', ...
                        c_bridge', max_c_nextperiod_nobridge');
    tab_c_max.Properties.VariableNames = {'coh_today', 'aprime', 'aprime_c_cost_total', 'aprime_nobridge',...
                        'c_bridge', 'max_c_nextperiod_nobridge'};
    disp(tab_c_max);    
end


%% Process Results
fl_for_borr = zeros(size(c_infonly));
fl_for_save = zeros(size(c_infonly));
fl_inf_borr = zeros(size(c_infonly));
fl_coh_add = zeros(size(c_infonly));

idx_infonly = (c_infonly==max_c_nextperiod_nobridge);
idx_infforb = (c_infforb==max_c_nextperiod_nobridge);
idx_forbrsv = (c_forbrsv==max_c_nextperiod_nobridge);

fl_inf_borr(idx_infonly) = b_infonly_inf(idx_infonly);
fl_coh_add(idx_infonly) = c_infonly(idx_infonly);

fl_for_borr(idx_infforb) = b_infforb_for(idx_infforb);
fl_inf_borr(idx_infforb) = b_infforb_inf(idx_infforb);
fl_coh_add(idx_infforb) = c_infforb(idx_infforb);

fl_for_borr(idx_forbrsv) = b_forbrsv_brr(idx_forbrsv);
fl_for_save(idx_forbrsv) = b_forbrsv_sav(idx_forbrsv);
fl_coh_add(idx_forbrsv) = c_forbrsv(idx_forbrsv);

end
