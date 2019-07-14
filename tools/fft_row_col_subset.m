%% Given Row and Column Counts, Get Subset of Rows and Columns for Display
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [ar_it_cols, ar_it_rows] = fft_row_col_subset(varargin)
%% FFT_ROW_COL_SUBSET
%
% @example
%
%    [it_row_n, it_col_n] = size(mt_z_trans);
%    [ar_it_cols, ar_it_rows] = fft_row_col_subset(it_col_n, it_col_n_keep, it_row_n, it_row_n_keep);    
%    cl_st_full_rowscols = cellstr([num2str(ar_z_r_borr_mesh_wage', 'zr=%3.2f;'), ...
%                                   num2str(ar_z_wage_mesh_r_borr', 'zw=%3.2f')]);
%    tb_mt_z_trans = array2table(round(mt_z_trans(ar_it_rows, ar_it_cols),6));
%    cl_col_names = strcat('zi=', num2str(ar_it_cols'), ':', cl_st_full_rowscols(ar_it_cols));
%    cl_row_names = strcat('zi=', num2str(ar_it_rows'), ':', cl_st_full_rowscols(ar_it_cols));
%    tb_mt_z_trans.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
%    tb_mt_z_trans.Properties.RowNames = matlab.lang.makeValidName(cl_row_names);
%    
%    disp(size(tb_mt_z_trans));
%    disp(tb_mt_z_trans);
%

%% Default
% use binomial as test case, z maps to binomial win prob, remember binom
% approximates normal.

if (isempty(varargin))
    
    clear all;
    close all;
    
    [it_col_n, it_col_n_keep, it_row_n, it_row_n_keep] ...
        = deal(100, 4, 37, 4);
    bl_display_addzero = true;
    
else
    
    [it_col_n, it_col_n_keep, it_row_n, it_row_n_keep] = varargin{:};
    bl_display_addzero = false;
    
end

%% Add in Zero

if (it_col_n >= it_col_n_keep)
    
    it_lower_end = floor(it_col_n_keep/2);
    it_upper_str = ((it_col_n)-floor(it_col_n_keep/2)+1);
        
    ar_it_cols = (1:1:it_lower_end);
    ar_it_cols = [ar_it_cols it_upper_str:1:(it_col_n)];

    if (rem(it_col_n_keep, 2) ~= 0)
        % add mid point in 
        ar_it_cols = [ar_it_cols round(it_col_n/2)];
    end
    
else
    ar_it_cols = 1:1:it_col_n;
end
ar_it_cols = sort(unique(ar_it_cols));

if (it_row_n >= it_row_n_keep)
    
    it_lower_end = floor(it_row_n_keep/2);
    it_upper_str = ((it_row_n)-floor(it_row_n_keep/2)+1);
    
    ar_it_rows = (1:1:it_lower_end);
    ar_it_rows = [ar_it_rows it_upper_str:1:(it_row_n)];

    if (rem(it_row_n_keep, 2) ~= 0)
        % add mid point in 
        ar_it_rows = [ar_it_rows round(it_row_n/2)];
    end
    
else
    ar_it_rows = 1:1:it_row_n;
end
ar_it_rows = sort(unique(ar_it_rows));


%% Display
if (bl_display_addzero)
    
    disp('ar_it_cols');
    disp(ar_it_cols);
    
    disp('ar_it_rows_iter');
    disp(ar_it_rows);
    
end
end
