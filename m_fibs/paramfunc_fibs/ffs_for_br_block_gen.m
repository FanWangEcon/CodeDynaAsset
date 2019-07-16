%% Generate Formal Borrowing Menu
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [ar_forbrblk, ar_forbrblk_r] = ffs_for_br_block_gen(varargin)
%% FFS_FOR_BR_BLOCK_GEN formal borrowing blocks
% Grid of formal borrowing points with associated interest rates. Policies
% could shift this. This function could be invoked to generate multiple
% formal borrowing grids from multiple formal borrowing providers, and the
% overall formal borrowing choice set is a combination of these grids. 
%
% @param fl_r_fbr boolean single formal borrowing interest rate
%
% @param st_forbrblk_type string different formal borrowing grid
% structures, some for testing/illustration purposes. Estimation uses real
% observed blocks from data.
%
% @param st_forbrblk_type string different formal borrowing grid
% structures, some for testing/illustration purposes. Estimation uses real
% observed blocks from data.
%
% # *unif*: uniform grid equal gap from min to
% # *seg3*: three segment blocks, the first 1/3 gap of the third segement,
% respect min and max from fl_forbrblk_brmost and fl_forbrblk_brleast,
% smallest borrowing segment gap is fl_forbrblk_gap
%
% @param fl_forbrblk_brmost float maximum formal borrowing allowed on grid,
% , might not not be relevant for all _st_forbrblk_type_
%
% @param fl_forbrblk_brleast float minimum formal borrowing allowed on
% grid, might not be relevant for all _st_forbrblk_type_
%
% @param fl_forbrblk_gap float formal borrowing grid gap, means different
% things for different _st_forbrblk_type_
%
% @return ar_forbrblk array array of formal borrowing grid points. should
% be ordered from lowest to highest including zero at the top. i.e.: [-4,
% -2, -1, 0]. These are assumed to be in principles units always.
% PRINCIPLES only, do not include interest rates.
%
% @return ar_forbrblk_r array interest rates associated with equal-length
% _ar_forbrblk_. ar_forbrblk_r multiplied element-wise with ar_forbrblk
% will provide the principle + interest rates associated with each formal
% borrowing grid elements. 
%
% @example
%
%   [ar_forbrblk, ar_forbrblk_r] = ...
%        ffs_for_br_block_gen(fl_r_fbr, st_forbrblk_type, fl_forbrblk_brmost, fl_forbrblk_brleast, fl_forbrblk_gap);
%
% @seealso
%
% * Formal Borrowing Grid: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_for_br_block_gen.html ffs_for_br_block_gen>
% * Informal Interest Rates: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_r_inf.html ffs_r_inf>
% * Match Borrowing to Formal Grid: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_for_br_block_match.html ffs_for_br_block_match>
% * Optimize Formal and Informal, Borrowing and Savings Joint Choices: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_min_c_cost.html ffs_fibs_min_c_cost>
% * Bridge Loan: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_inf_bridge.html ffs_fibs_inf_bridge>
% * Overall Optimization: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_min_c_cost_bridge.html ffs_fibs_min_c_cost_bridge>
% * Discrete Choices: <https://fanwangecon.github.io/CodeDynaAsset/m_fibs/paramfunc_fibs/html/ffs_fibs_identify_discrete.html ffs_fibs_identify_discrete>
%


%% Default

fl_r_fbr = 0.045;
[st_forbrblk_type, fl_forbrblk_brmost, fl_forbrblk_brleast, fl_forbrblk_gap] = deal('unif', -19, -1, -0.5);
bl_display_forbrblock = false;
default_params = {fl_r_fbr st_forbrblk_type fl_forbrblk_brmost fl_forbrblk_brleast fl_forbrblk_gap ...
    bl_display_forbrblock};

%% Parse Parameters

% numvarargs is the number of varagin inputted
[default_params{1:length(varargin)}] = varargin{:};
[fl_r_fbr, st_forbrblk_type, fl_forbrblk_brmost, fl_forbrblk_brleast, fl_forbrblk_gap, ...
    bl_display_forbrblock] = default_params{:};

%% Set Formal Borrowing Grid Uniform

if (strcmp(st_forbrblk_type, 'unif'))
    ar_forbrblk = fl_forbrblk_brleast:fl_forbrblk_gap:fl_forbrblk_brmost;
end

%% Set Formal Borrowing Grid Three Segment
% An example grid with three segments, the first segement has
% fl_forbrblk_gap as gap, the second segment gap is twice, the third
% segement the gap is 3 times as large as fl_forbrblk_gap. Respect
% fl_forbrblk_brmost fl_forbrblk_brleast, which means divided into three
% parts based on them as much as possible

if (strcmp(st_forbrblk_type, 'seg3'))
    fl_most_least_gap = (fl_forbrblk_brmost - fl_forbrblk_brleast);
    fl_most_least_seg3_interval = fl_most_least_gap/3;
    
    ar_seg1 = fl_forbrblk_brleast:fl_forbrblk_gap:fl_most_least_seg3_interval;
    ar_seg2 = max(ar_seg1):(fl_forbrblk_gap*2):fl_most_least_seg3_interval*2;
    ar_seg3 = max(ar_seg2):(fl_forbrblk_gap*3):fl_forbrblk_brmost;
    
    ar_forbrblk =[ar_seg1, ar_seg2, ar_seg3];
end

%% Set Formal Borrowing Grid Three Segment
% An example grid with three segments, the first segement has
% fl_forbrblk_gap as gap, the second segment gap is twice, the third
% segement the gap is 3 times as large as fl_forbrblk_gap. Respect
% fl_forbrblk_brmost fl_forbrblk_brleast, which means divided into three
% parts based on them as much as possible

if (strcmp(st_forbrblk_type, 'testfx'))
    fl_most_least_gap = (fl_forbrblk_brmost - fl_forbrblk_brleast);
    fl_most_least_seg3_interval = fl_most_least_gap/3;
    
    ar_seg1 = fl_forbrblk_brleast:fl_forbrblk_gap:fl_most_least_seg3_interval;
    ar_seg2 = max(ar_seg1):(fl_forbrblk_gap*2):fl_most_least_seg3_interval*2;
    ar_seg3 = max(ar_seg2):(fl_forbrblk_gap*3):fl_forbrblk_brmost;
    
    ar_forbrblk =[-3 -6];
end


%% Sort Borrowing Blocks

% add zero because formal borrowing = 0 should be an option
ar_forbrblk = [-0 ar_forbrblk];
ar_forbrblk = sort(unique([ar_forbrblk]));

%% Generate Corresponding Interest Rates Arrays
ar_forbrblk_r = fl_r_fbr + zeros(size(ar_forbrblk));

if (bl_display_forbrblock)
    disp('ar_forbrblk');
    disp(ar_forbrblk);
    disp('ar_forbrblk_r');
    disp(ar_forbrblk_r);
end

end
