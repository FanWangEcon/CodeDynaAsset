%%
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

function [ar_a_grid_ceil_principle, ar_a_grid_ceil_wthr, ...
    ar_a_grid_floor_principle, ar_a_grid_floor_wthr] = ffs_for_br_block_match(varargin)
%% FFS_FOR_BR_BLOCK_MATCH formal borrowing blocks
% Find Value just below or above each element of *ar_a* from *ar_forbrblk*.
% a vector of grid points, find for each element of ar_a the element of
% ar_forbrblk that is just above or just below.
%
% # _ar_a_: $a'_i$ where $i \in (1,...,N)$
% # _ar_forbrblk_: $grid_j$ where $j\in (1,...,M)$
%
% $$ForCEIL_i = argmin_{j}(ForGrid_j - a'_i | grid_j > a'_i)$$
%
% $ForCEIL_i$ is the level of formal borrowing if joint formal + informal
% borrowing is chosen.
%
% $$ForFLOOR_i = argmax_{j}(ForGrid_j - a'_i | grid_j <= a'_i)$$
%
% $ForFLOOR_i$ is the level of formal borrowing if joint formal borrowing +
% informal savings is chosen.
%
% @param ar_a boolean N by 1 single formal borrowing interest rate
%
% @param ar_forbrblk 1 by M array array of formal borrowing grid points
%
% @param ar_forbrblk_r array interest rates associated with equal-length
% _ar_forbrblk_.
%
% @param bl_b_is_principle boolean solving with aggregate savings as
% savings + debt principles + interests, or just principles no interests.
% if true, principels only, no interests.
%
% @return ar_a_grid_ceil_principle array N by 1 Solution to:
%
% $$min_{j}(ForGrid_j - a'_i | grid_j > a'_i)$$
%
% @return ar_a_grid_ceil_wthr array _ar_a_grid_ceil_principle_ with
% interest rates specified to each borrowing formal level added
%
% @return ar_a_grid_floor_principle array N by 1 element of the *ar_forbrblk* vector that are
% the elements right above each eelemnt of ar_a. Solution to:
%
% $$max_{j}(ForGrid_j - a'_i | grid_j <= a'_i)$$
%
% @return ar_a_grid_floor_principle array _ar_a_grid_floor_principle_ with
% interest rates specified to each borrowing formal level added
%
% @example
%
%   [ar_a_grid_ceil, ar_a_grid_floor] = ...
%        ffs_for_br_block_match(ar_a, ar_forbrblk, ar_forbrblk_r, bl_b_is_principle);
%


%% Default

% array of a choices
% ar_a could be principles + interests, or principles only
ar_a = -sort(rand([10,1])*20);

% use defaults from block gen
[ar_forbrblk, ar_forbrblk_r] = ffs_for_br_block_gen();

% if bl_b_is_principle is true, b is principles only, no interests.
% bl_b_is_principle = false is the case for models like *abz* without
% interpolation over cash-on-hand
bl_b_is_principle = false;

% Display
bl_display_brblockmatch = false;

default_params = {ar_a ar_forbrblk ar_forbrblk_r bl_b_is_principle bl_display_brblockmatch};

%% Parse Parameters

% numvarargs is the number of varagin inputted
[default_params{1:length(varargin)}] = varargin{:};
[ar_a, ar_forbrblk, ar_forbrblk_r, bl_b_is_principle, bl_display_brblockmatch] = default_params{:};

%% Adjust Inputs t
% if bl_b_is_principle, then principle, with the assumption that
% ar_forbrblk. If bl_b_is_principle is false, that means the ar_a vector is
% principle and interest rates. Hence, need to convert ar_forbrblk which
% are principles to interests plus principles to be on the same scale as
% ar_a.

if (~bl_b_is_principle)
    ar_forbrblk_use = ar_forbrblk.*(1+ar_forbrblk_r);
else
    ar_forbrblk_use = ar_forbrblk;
end

%% Show Details Step by Step
if (bl_display_brblockmatch)

    % show borrowing array
    disp('ar_a')
    disp(ar_a)

    % show borrowing formal blocks/grids
    disp('ar_forbrblk_use and ar_forbrblk')
    disp([ar_forbrblk_use;ar_forbrblk]')

    % all combination division
    disp('mt_a_dvd_grid = (ar_a./ar_forbrblk_use)')
    mt_a_dvd_grid = (ar_a./ar_forbrblk_use)

    % ceiling for each
    disp('(mt_a_dvd_grid >= 1)')
    (mt_a_dvd_grid >= 1)

    % If ceiling exists and cloest ceiling index
    % min_{j}( ar_forbrblk[j] - ar_a[i] | ar_forbrblk[j] > ar_a[i])
    disp('[~, ar_max_a_on_grid_idx] = max((mt_a_dvd_grid >= 1),[], 2)')
    [~, ar_max_a_on_grid_idx] = max((mt_a_dvd_grid >= 1),[], 2)

    % ar_forbrblk[argmin_{j}( ar_forbrblk[j] - ar_a[i] | ar_forbrblk[j] > ar_a[i])]
    disp('ar_a_grid_ceil = ar_forbrblk_use(ar_max_a_on_grid_idx)')
    ar_a_grid_ceil = ar_forbrblk_use(ar_max_a_on_grid_idx)
    % ar_a_grid_ceil(ar_max_a_on_grid_idx == 1) = ar_forbrblk(0)

    % now floor, just one index less
    disp('ar_a_grid_floor = ar_forbrblk_use(max(ar_max_a_on_grid_idx - 1, 1))')
    ar_a_grid_floor = ar_forbrblk_use(max(ar_max_a_on_grid_idx - 1, 1))
    % ar_a_grid_floor(ar_max_a_on_grid_idx == 1) =

    % Dispaly
    tab_matched_grid = table(ar_a, ar_a_grid_floor', ar_a_grid_ceil');
    tab_matched_grid.Properties.VariableNames = {'ar_a','ar_a_grid_floor','ar_a_grid_ceil'};
    disp('ar_a_grid_floor: for borrow + for save')
    disp('ar_a_grid_ceil: for + inf borrow')
    disp(tab_matched_grid);
end

%% Standard Quicker Solve

% Get Index
[~, ar_max_a_on_grid_idx] = max(((ar_a./ar_forbrblk_use) >= 1),[], 2);

% Get Values
if (~bl_b_is_principle)

    % Borrowing borrowing points, following formal grids, but add interests
    ar_a_grid_ceil_wthr = ar_forbrblk_use(ar_max_a_on_grid_idx)';
    ar_a_grid_floor_wthr = ar_forbrblk_use(max(ar_max_a_on_grid_idx - 1, 1))';

    % Principles only
    ar_a_grid_ceil_principle = ar_forbrblk(ar_max_a_on_grid_idx)';
    ar_a_grid_floor_principle = ar_forbrblk(max(ar_max_a_on_grid_idx - 1, 1))';

else

    % Borrowing borrowing points, following formal grids, but add interests
    ar_a_grid_ceil_wthr = ...
        (ar_forbrblk_use(ar_max_a_on_grid_idx).*ar_forbrblk_r(ar_max_a_on_grid_idx))';
    ar_a_grid_floor_wthr = ...
        (ar_forbrblk_use(max(ar_max_a_on_grid_idx - 1, 1)).*ar_forbrblk_r(max(ar_max_a_on_grid_idx - 1, 1)))';

    % Principles only, note ar_forbrblk_use = ar_forbrblk
    ar_a_grid_ceil_principle = ar_forbrblk_use(ar_max_a_on_grid_idx)';
    ar_a_grid_floor_principle = ar_forbrblk_use(max(ar_max_a_on_grid_idx - 1, 1))';

end
