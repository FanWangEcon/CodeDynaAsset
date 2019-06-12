function [fl_for_borr, fl_for_save, fl_inf_borr, fl_coh_add] = ...
    borrtranslate_ar(fl_r_inf, fl_r_fsv, fl_r_fbr, fl_for_br_block, aprime)

% %% Parameters
% % The informal interest rate
% fl_r_inf = 0.10;
%
% % The formal interest rate save for borrowing/savings
% fl_r_fsv = 0.01;
% fl_r_fbr = 0.06;
%
% % current aprime amount = principles + interest
% aprime = -5;
%
% % Formal borrowing blocks, minimal formal borrowing sizes
% fl_for_br_block = -1;

%% Consumption Values
% Informal Borrowing only
b_infonly_inf = (aprime./(1+fl_r_inf));
c_infonly = (-1).*b_infonly_inf;

% Consumption Values Formal + Informal Borrowing Jointly
for_br_block_dvd = (aprime./(1+fl_r_fbr))./fl_for_br_block;
for_br_block_dvd_floor = abs(floor(for_br_block_dvd).*fl_for_br_block);
for_br_block_dvd_ceil = abs(ceil(for_br_block_dvd).*fl_for_br_block);

b_infforb_inf = ((aprime + (1+fl_r_fbr).*for_br_block_dvd_floor)./(1+fl_r_inf));
b_infforb_for = -(for_br_block_dvd_floor);
c_infforb =  -1*(b_infforb_inf + b_infforb_for);

% Consumption Values Formal Borrowing + Formal Savings
b_forbrsv_sav = ((aprime + (1+fl_r_fbr).*for_br_block_dvd_ceil)./(1+fl_r_fsv));
b_forbrsv_brr = -(for_br_block_dvd_ceil);
c_forbrsv =  -1*(b_forbrsv_sav + b_forbrsv_brr);

%% Processing
fl_for_borr = zeros(size(c_infonly));
fl_for_save = zeros(size(c_infonly));
fl_inf_borr = zeros(size(c_infonly));
fl_coh_add = zeros(size(c_infonly));

% Trying to minimize consumption costs
max_c = max([c_infonly; c_infforb; c_forbrsv],[], 1);

idx_infonly = (c_infonly==max_c);
idx_infforb = (c_infforb==max_c);
idx_forbrsv = (c_forbrsv==max_c);

fl_inf_borr(idx_infonly) = b_infonly_inf(idx_infonly);
fl_coh_add(idx_infonly) = c_infonly(idx_infonly);

fl_for_borr(idx_infforb) = b_infforb_for(idx_infforb);
fl_inf_borr(idx_infforb) = b_infforb_inf(idx_infforb);
fl_coh_add(idx_infforb) = c_infforb(idx_infforb);

fl_for_borr(idx_forbrsv) = b_forbrsv_brr(idx_forbrsv);
fl_for_save(idx_forbrsv) = b_forbrsv_sav(idx_forbrsv);
fl_coh_add(idx_forbrsv) = c_forbrsv(idx_forbrsv);


end
