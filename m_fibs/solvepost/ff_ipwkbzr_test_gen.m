%% Test Risky + Safe Asset (Save + Borr + FIBS + RShock) Interp-Perc Main Function
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function tb_outcomes = ff_ipwkbzr_test_gen(varargin)
%% FF_IPWKBZR_TESTMAIN main formal informal borrow save tester
% Testing Program

%% Default Parameter Loops

% Loops to Vary for Simulations Below
it_size_type = 1;
ar_it_test_grp = [3, 8, 9];
it_simu_vec_len= 3;
bl_simu_cross = true; % if false, simu full grid

%% Default Parameters

% Base parameters call
it_param_set = 9;
[param_map, support_map] = ffs_ipwkbzr_fibs_set_default_param(it_param_set);
support_map('bl_timer') = true;
support_map('bl_mat_test') = true;
st_matimg_path_root = support_map('st_matimg_path_root');
% test_borinf for default: ar_it_test_grp = [3, 8, 9]
support_map('st_mat_test_path') = [st_matimg_path_root '/test/ff_ipwkbzr_ds_vecsv/test_forinf/mat/'];

%% Array Parameters

ls_st_param_key = {'fl_b_bd', 'fl_c_min', 'fl_z_r_infbr_poiss_mean', ...
    'fl_beta', 'fl_crra', 'fl_z_wage_rho', 'fl_z_wage_sig', ...
    'fl_r_fbr', 'fl_forbrblk_gap'};

ls_st_param_desc = {'Borrow Bound', 'Minimum Consumption', 'Mean Inf Interest', ...
    'Discount', 'CRRA', 'Shock Persistence', 'Shock SD', ...
    'Formal Borrow R', 'For Borrow Menu Gap'};

param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
param_tstar_map('it_z_wage_n') = [5];

% Length all reset later based on the min and max here
param_tstar_map(ls_st_param_key{1}) = linspace(-20, -5, it_simu_vec_len);
param_tstar_map(ls_st_param_key{2}) = linspace(0.1, 0.001, it_simu_vec_len);
param_tstar_map(ls_st_param_key{4}) = linspace(0.94, 0.98, it_simu_vec_len);
param_tstar_map(ls_st_param_key{5}) = linspace(1, 2, it_simu_vec_len);
param_tstar_map(ls_st_param_key{6}) = linspace(0.65, 0.95, it_simu_vec_len);
param_tstar_map(ls_st_param_key{7}) = linspace(0.05, 0.35, it_simu_vec_len);

param_tstar_map(ls_st_param_key{3}) = linspace(1, 20, it_simu_vec_len);
param_tstar_map(ls_st_param_key{8}) = linspace(0.030, 0.080, it_simu_vec_len);
param_tstar_map(ls_st_param_key{9}) = linspace(-0.1, -2.5, it_simu_vec_len);


%% Default
default_params = {it_size_type ar_it_test_grp ...
    param_map support_map param_tstar_map};

%% Parse Inputs
% override default set above if any parameters are updated

params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};
it_size_type = default_params{1};
ar_it_test_grp = default_params{2};
param_map = [param_map; default_params{3}];
support_map = [support_map; default_params{4}];
param_tstar_map = [param_tstar_map; default_params{5}];

% Generate additional Elements of param_map for container map display
param_map('it_size_type') = it_size_type;
param_map('ar_it_test_grp') = ar_it_test_grp;
ls_st_param_key_subset = ls_st_param_key(ar_it_test_grp);
param_map('ls_st_param_key_subset') = ls_st_param_key_subset;

%% Store Mat Strings
% support_map

support_map('st_mat_test_prefix') = [''];
support_map('st_mat_test_name_main') = ['res'];
support_map('st_mat_test_suffix') = ['g' strrep(num2str(ar_it_test_grp), '  ', '') ...
    '_c' num2str(bl_simu_cross) '_t' num2str(it_size_type)];

%% Set Solve Sizes

if (it_size_type == 1)
    
    % Basic Test Run
    param_map('it_w_perc_n') = 10;
    param_map('it_ak_perc_n') = param_map('it_w_perc_n');
    param_map('it_coh_bridge_perc_n') = param_map('it_w_perc_n');
    
    param_map('fl_coh_interp_grid_gap') = 2;
    param_map('it_c_interp_grid_gap') = 10^-1;
    param_map('fl_w_interp_grid_gap') = 2;
    
    param_map('it_z_wage_n') = 2;
    param_map('fl_z_r_infbr_n') = 2;
    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');
    
    param_map('it_maxiter_val') = 5;
    
elseif (it_size_type == 2)
    
    % Usable Run
    param_map('it_maxiter_val') = 50;
    
elseif (it_size_type == 3)
    
    % Full Run
    
elseif (it_size_type == 4)
    
    % Denser Run
    param_map('it_w_perc_n') = 100;
    param_map('it_ak_perc_n') = param_map('it_w_perc_n');
    % ususally divide by 5
    param_map('it_coh_bridge_perc_n') = round(param_map('it_w_perc_n')/3);
    
end

%% Parase Preference and Shock Parameters

% support map
params_group = values(support_map, {'bl_mat_test', 'st_mat_test_path', 'st_mat_test_prefix', 'st_mat_test_name_main', 'st_mat_test_suffix'});
[bl_mat_test, st_mat_test_path, st_mat_test_prefix, st_mat_test_name_main, st_mat_test_suffix] = params_group{:};

%% Display support_map

fft_container_map_display(support_map);
fft_container_map_display(param_map);
fft_container_map_display(param_tstar_map);

%% Initialize Storage

disp('---------------------------');
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp(['Vary These Parameters:' ls_st_param_key_subset]);
disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
disp('---------------------------');

% cl_tb_outcomes_meansdperc_wthinfo = cell([length(ar_it_test_grp), 1]);
tb_outcomes = [];

%% Simulate 1: Cross Simulation, fix one point, extend parameters in each direction
% Given X and Y parameters, simulate along an array of x values fixing y
% value, then simulate along an array of y values fixing x value. Along the
% x and y grid, this produces a _cross_ of simulated points. The
% intersection of the array arrays, the center of the cross could be some
% benchmark parameter value, or perhaps some estimated/calibrated parameter
% value point.
%

if (bl_simu_cross)
    for it_pcombi_ctr = 1:length(ar_it_test_grp)
        
        it_test_grp = ar_it_test_grp(it_pcombi_ctr);
        
        % Display Current Parameter been Varied        
        % Parameter Key for Paraemter Getting Updated
        st_param_key = ls_st_param_key{it_test_grp};
        st_param_desc = ls_st_param_desc{it_test_grp};
        ar_param_values = param_tstar_map(st_param_key);
        fl_param_val_benchmark = param_map(st_param_key);
        
        disp('---------------------------');
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        disp(['Vary ' st_param_desc ' (' st_param_key '): ' num2str(ar_param_values)]);
        for it_param_all_ctr = 1:length(ls_st_param_key)
            disp([ ls_st_param_key{it_param_all_ctr} ':' num2str(param_map(ls_st_param_key{it_param_all_ctr}))]);
        end
        disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
        
        % Simulate Model over Parameter Array Values        
        for it_cur_param = 1:1:length(ar_param_values)
            
            % Adjust Value for Current Parameter been Varied
            fl_param_val = ar_param_values(it_cur_param);
            param_map(st_param_key) = fl_param_val;
            disp(['xxxxx ' st_param_key ' = ' num2str(fl_param_val) ' xxxxx']);
            
            % Simulate Model
            ar_simu_info = [it_test_grp, it_cur_param, cell2mat(values(param_map, ls_st_param_key_subset))];           
            cl_col_names = [{'it_test_grp', 'it_cur_param'} ls_st_param_key_subset];
            [tb_outcomes_simu, variablenames] = simu_model_gen_stats(param_map, support_map, ar_simu_info, cl_col_names);
            tb_outcomes_simu.Properties.RowNames = ...
                strcat(variablenames, '_p', num2str(it_pcombi_ctr), 'v', num2str(it_cur_param));
            tb_outcomes_simu = addvars(tb_outcomes_simu, variablenames, 'Before', 1);
            
            % Combine Results from Different Parameters in Shared Table
            if (it_pcombi_ctr == 1 && it_cur_param == 1)
                tb_outcomes = tb_outcomes_simu;
            else
                tb_outcomes = [tb_outcomes; tb_outcomes_simu];
            end
            
        end        
        % Reset Base Parameters, parameters already grabbed out, updating
        param_map(st_param_key) = fl_param_val_benchmark;        
    end
end

%% Simulate 2: Full Grid Simulation, Simulate Along Full Grid
% To explore the effects of parameters on model outcomes, simulate the
% model along full grids. Given X and Y parameters, this means simulate at
% all possible combinations of X and Y arrays.
%

% Get Arrays to be Meshed in Cells
cl_ar_param_subset_values = values(param_tstar_map, ls_st_param_key_subset);

% Generate all possible combinations of parameters subsets
cl_mt_all = cl_ar_param_subset_values;
[cl_mt_all{:}] = ndgrid(cl_ar_param_subset_values{:});
mt_param_vals_combi = cell2mat(cellfun(@(m) m(:), cl_mt_all, 'uni', 0));

% Sizes
it_test_combi_n = size(mt_param_vals_combi,1);

% Convert from Matrix to Table for Clarity
tb_pvals_combi = array2table(mt_param_vals_combi);
tb_pvals_combi.Properties.VariableNames = ls_st_param_key(ar_it_test_grp);
tb_pvals_combi.Properties.RowNames = strcat('j=', string(1:size(mt_param_vals_combi,1)));
clear mt_param_vals_combi;

% Display
disp(head(tb_pvals_combi, 10));
disp(tail(tb_pvals_combi, 10));

if (~bl_simu_cross)
    for it_pcombi_ctr = 1:it_test_combi_n
        
        tb_row_param_adj_cur = tb_pvals_combi(it_pcombi_ctr, :);
        
        disp(['xxxxx Shift to: xxxxx']);
        disp(tb_row_param_adj_cur)
        
        % Display Current Parameter been Varied        
        for st_param_key = ls_st_param_key_subset
            param_map(st_param_key{1}) = tb_row_param_adj_cur{1, st_param_key};
        end
            
        % Simulate Model        
        ar_simu_info = [it_pcombi_ctr, cell2mat(values(param_map, ls_st_param_key_subset))];           
        cl_col_names = [{'it_pcombi_ctr'} ls_st_param_key_subset];        
        [tb_outcomes_simu, variablenames] = simu_model_gen_stats(param_map, support_map, ar_simu_info, cl_col_names);
        tb_outcomes_simu.Properties.RowNames = strcat(variablenames, '_v', num2str(it_pcombi_ctr));
        tb_outcomes_simu = addvars(tb_outcomes_simu, variablenames, 'Before', 1);

        % Combine Results from Different Parameters in Shared Table
        if (it_pcombi_ctr == 1)
            tb_outcomes = tb_outcomes_simu;
        else
            tb_outcomes = [tb_outcomes; tb_outcomes_simu];
        end

    end    
end

%% Save Mat
if (bl_mat_test)
    clear armt_map result_map
    if ~exist(st_mat_test_path,'dir'); mkdir(st_mat_test_path); end
    st_file_name = [st_mat_test_prefix st_mat_test_name_main st_mat_test_suffix];
    save(strcat(st_mat_test_path, st_file_name));
end

end

%% Simulate Model given Parameter and Process Results
function [tb_outcomes_simu, variablenames] = simu_model_gen_stats(param_map, support_map, ar_simu_info, cl_col_names)

    ls_st_param_key_subset = param_map('ls_st_param_key_subset');
    
    %% Reset Parameters that are determined by other parameters
    % in case other parameters have been changed
    param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_infbr_n');
    it_w_perc_n = param_map('it_w_perc_n');
    it_z_n = param_map('it_z_n');
    disp(['xxxxx it_w_perc_n = ' num2str(it_w_perc_n) ', it_z_n = ' num2str(it_z_n) ' xxxxx']);

    %% Simulate Model

    [armt_map, func_map] = ffs_ipwkbzr_fibs_get_funcgrid(param_map, support_map);
    result_map = ff_ipwkbzr_fibs_vf_vecsv(param_map, support_map, armt_map, func_map);
    result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map);

    %% Store Results to Table

    % Append Simulation Loop Statistics
    tb_outcomes_meansdperc = result_map('tb_outcomes_meansdperc');
    tb_simu_info = array2table(ar_simu_info + zeros([size(tb_outcomes_meansdperc,1), length(ar_simu_info)]));    
    tb_simu_info.Properties.VariableNames = cl_col_names;
    variablenames = tb_outcomes_meansdperc.Properties.RowNames;
    tb_simu_info.Properties.RowNames = variablenames;
    
    % Table
    tb_outcomes_simu = [tb_simu_info tb_outcomes_meansdperc];
end