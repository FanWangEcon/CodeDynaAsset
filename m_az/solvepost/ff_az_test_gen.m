%% Simulate Savings Model Along Various Parameter Arrays
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*
  
%%
function [tb_outcomes, support_map, param_desc_map] = ff_az_test_gen(varargin)
%% FF_AZ_TEST_GEN post solution simulation
% Simulate the model along various dimensions
%
% @param bl_simu_cross boolean cross vs gridd simulation. cross: with
% (x,y), vary X fixing y, vary Y fixing x. grid: all (x,y) \in (X,Y)
%
% @param it_size_type integer: 
%  param_map support_map param_tstar_map param_desc_map
% # it_size_type = 1 is quick grid
% # it_size_type = 2 is standard grid
% # it_size_type = 3 is denser grid
%
% @param cl_st_param_keys cell string cell array container parameters that
% we are simulating over
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param param_tstar_map container map of arrays with keys that are
% parameter keys. This could be specified outside with array values to
% override defaults here. 
%
% @param param_desc_map container map of strings for each parameter key.
%
% @param tb_outcomes table table of simulation outcomes for various
% statistics for each set of parameters. Columns are statistics, rows are
% outcome variables, groups of rows are different simulations.
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html ff_az_vf_vecsv>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vecsv.html ff_az_ds_vecsv>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html ffs_az_set_default_param>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_get_funcgrid.html ffs_az_get_funcgrid>
%
% @seealso
%
% * _SPEED_ savings only overall benchmark speed testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_speed/html/fsi_az_ds_vecsv_speed.html fsi_az_ds_vecsv_speed>
% * _PREFERENCE_ savings only preference testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_pref/html/fsi_az_ds_vecsv_pref.html fsi_az_ds_vecsv_pref>
% * _PREFERENCE_ savings only preference testing cross:
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_pref/html/fsi_az_ds_vecsv_pref_cross.html fsi_az_ds_vecsv_pref_cross>
% * _SHOCK_ savings only shock testing: <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_shock/html/fsi_az_ds_vecsv_shock.html fsi_az_ds_vecsv_shock>
% * _SHOCK_ savings only shock testing cross:
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_shock/html/fsi_az_ds_vecsv_shock_cross.html fsi_az_ds_vecsv_shock_cross>
% * _PRICE_ savings only wage and interest rate testing cross: adjust wage and savings rate
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/test/ff_az_ds_vecsv/test_price/html/fsi_az_ds_vecsv_price_cross.html fsi_az_ds_vecsv_price_cross>
%

%% Default Parameter Loops

% Loops to Vary for Simulations Below
it_size_type = 2;
cl_st_param_keys = {'fl_crra', 'fl_beta'};
bl_simu_cross = true; % if false, simu full grid

%% Default Parameters

% Base parameters call
it_param_set = 9;
[param_map, support_map] = ffs_az_set_default_param(it_param_set);

% Timer
support_map('bl_timer') = true;

% Mat Storage
support_map('bl_mat_test') = true;
st_matimg_path_root = support_map('st_matimg_path_root');
% test_borinf for default: cl_st_param_keys = [3, 8, 9]
support_map('st_mat_test_path') = [st_matimg_path_root '/test/ff_az_ds_vecsv/mat/'];

% Generate or Load
% if bl_replacefile is true, that means even if file already exists,
% generaet new. bl_replacefile is false, then if false does not exist
% generate, if file exists does not generate.
support_map('bl_replacefile') = true;
support_map('bl_display_simu_stats') = true;

%% Array Parameters

ls_st_param_key = {'fl_crra', 'fl_beta', ...
                    'fl_w', 'fl_r_save', ...
                    'fl_z_rho', 'fl_z_sig', ...
                    'fl_a_max', 'it_z_n', 'it_a_n'};

param_tstar_map = containers.Map('KeyType','char', 'ValueType','any');
it_simu_vec_len = 5;
param_tstar_map(ls_st_param_key{1}) = linspace(1, 5, it_simu_vec_len);
param_tstar_map(ls_st_param_key{2}) = linspace(0.87, 0.97, it_simu_vec_len);
param_tstar_map(ls_st_param_key{3}) = linspace(1.1, 1.4, it_simu_vec_len);
param_tstar_map(ls_st_param_key{4}) = linspace(0.01, 0.04, it_simu_vec_len);
param_tstar_map(ls_st_param_key{5}) = linspace(0, 0.99, it_simu_vec_len);
param_tstar_map(ls_st_param_key{6}) = linspace(0.01, 0.5, it_simu_vec_len);
param_tstar_map(ls_st_param_key{7}) = linspace(50, 80, it_simu_vec_len);
param_tstar_map(ls_st_param_key{8}) = unique(round(linspace(5, 25, it_simu_vec_len)));
param_tstar_map(ls_st_param_key{9}) = unique(round(linspace(100, 2500, it_simu_vec_len)));

param_desc_map = containers.Map('KeyType','char', 'ValueType','any');
param_desc_map(ls_st_param_key{1}) = {'CRRA'};
param_desc_map(ls_st_param_key{2}) = {'Discount'};
param_desc_map(ls_st_param_key{3}) = {'Wage'};
param_desc_map(ls_st_param_key{4}) = {'Save Interest'};
param_desc_map(ls_st_param_key{5}) = {'Shock Persistence'};
param_desc_map(ls_st_param_key{6}) = {'Shock SD'};
param_desc_map(ls_st_param_key{7}) = {'Max Asset Bound'};
param_desc_map(ls_st_param_key{8}) = {'Shock Grid N'};
param_desc_map(ls_st_param_key{9}) = {'Asset Grid N'};

%% Default

default_params = {bl_simu_cross it_size_type cl_st_param_keys ...
                    param_map support_map param_tstar_map param_desc_map};

%% Parse Inputs
% override default set above if any parameters are updated

params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};

bl_simu_cross = default_params{1};
it_size_type = default_params{2};
cl_st_param_keys = default_params{3};

param_map = [param_map; default_params{4}];
support_map = [support_map; default_params{5}];
param_tstar_map = [param_tstar_map; default_params{6}];
param_desc_map = [param_desc_map; default_params{7}];

% Generate additional Elements of param_map for container map display
param_map('it_size_type') = it_size_type;
param_map('cl_st_param_keys') = cl_st_param_keys;
ar_param_keys_idx = cell2mat(cellfun(@(m) find(strcmp(ls_st_param_key, m)), ...
    cl_st_param_keys, 'UniformOutput', false));
param_map('ar_param_keys_idx') = ar_param_keys_idx;

%% Parse Generate or Get Current
% if bl_gen = true, generate, if bl_gen = false, get saved file and output.
% If bl_replacefile is true, that means the existing file must be replaced,
% which will set bl_gen to true

params_group = values(support_map, ...
    {'bl_replacefile', 'bl_display_simu_stats'});
[bl_replacefile, bl_display_simu_stats] = params_group{:};

%% Store Mat Strings
%
% # _it_total_length_: this is the total length of all selected arrays
% # _ar_param_keys_idx_: these are the linear index of the selected keys
% among all possible keys
%

it_total_length = sum(cell2mat(cellfun(@(m) length(param_tstar_map(m)), ...
                                cl_st_param_keys, 'UniformOutput', false)));

support_map('st_mat_test_prefix') = [''];
support_map('st_mat_test_name_main') = ['r'];
support_map('st_mat_test_suffix') = ['_g' strrep(num2str(ar_param_keys_idx), '  ', '') ...
                                     '_c' num2str(bl_simu_cross) 't' num2str(it_size_type) ...
                                     'l' num2str(it_total_length)];

%% Set Solve Sizes

if (it_size_type == 1)
    
    % Basic Test Run
    param_map('it_z_n') = 11;
    param_map('it_a_n') = 100;
        
elseif (it_size_type == 2)
    
    % Full Run
    
elseif (it_size_type == 3)
    
    % Denser Run
    param_map('it_z_n') = 27;
    param_map('it_a_n') = 2250;
    
end

%% Parase Preference and Shock Parameters

% support map
params_group = values(support_map, {'bl_mat_test', 'st_mat_test_path', 'st_mat_test_prefix', 'st_mat_test_name_main', 'st_mat_test_suffix'});
[bl_mat_test, st_mat_test_path, st_mat_test_prefix, st_mat_test_name_main, st_mat_test_suffix] = params_group{:};

%% Display support_map

fft_container_map_display(support_map);
fft_container_map_display(param_map);
fft_container_map_display(param_tstar_map);

%% Generate New or Obtain Existing Mat
% Mat file folder and mat file name

if ~exist(st_mat_test_path,'dir'); mkdir(st_mat_test_path); end
st_file_name = [st_mat_test_prefix st_mat_test_name_main st_mat_test_suffix '.mat'];
st_file_path_full = [st_mat_test_path, st_file_name];

bl_mat_exists = isfile(st_file_path_full);

if ( ~bl_replacefile && bl_mat_exists )
    
    st_loaded = load(st_file_path_full, 'tb_outcomes');
    tb_outcomes = st_loaded.tb_outcomes; 
    
else
    %% Initialize Storage
    
    disp('---------------------------');
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp(['Vary These Parameters:' cl_st_param_keys]);
    disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');
    disp('---------------------------');
    
    % cl_tb_outcomes_meansdperc_wthinfo = cell([length(cl_st_param_keys), 1]);
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
        
        for it_pcombi_ctr = 1:length(cl_st_param_keys)
            
            st_param_key = cl_st_param_keys{it_pcombi_ctr};
            
            % Display Current Parameter been Varied
            % Parameter Key for Paraemter Getting Updated
            ar_param_values = param_tstar_map(st_param_key);
            st_param_desc = param_desc_map(st_param_key);
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
                ar_simu_info = [find(strcmp(ls_st_param_key, st_param_key)), ...
                    it_cur_param, cell2mat(values(param_map, cl_st_param_keys))];
                cl_col_names = [{'it_test_grp', 'it_cur_param'} cl_st_param_keys];
                [tb_outcomes_simu, variablenames] = simu_model_gen_stats(param_map, support_map, ar_simu_info, cl_col_names);
                tb_outcomes_simu.Properties.RowNames = ...
                    strcat(variablenames, '_p', num2str(it_pcombi_ctr), 'v', num2str(it_cur_param));
                var_param_key = repmat({st_param_key}, [length(variablenames),1]);
                tb_outcomes_simu = addvars(tb_outcomes_simu, variablenames, var_param_key, 'Before', 1);
                
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
    
    if (~bl_simu_cross)
        % Get Arrays to be Meshed in Cells
        cl_ar_param_subset_values = values(param_tstar_map, cl_st_param_keys);
        
        % Generate all possible combinations of parameters subsets
        cl_mt_all = cl_ar_param_subset_values;
        [cl_mt_all{:}] = ndgrid(cl_ar_param_subset_values{:});
        mt_param_vals_combi = cell2mat(cellfun(@(m) m(:), cl_mt_all, 'uni', 0));
        
        % Sizes
        it_test_combi_n = size(mt_param_vals_combi,1);
        
        % Convert from Matrix to Table for Clarity
        tb_pvals_combi = array2table(mt_param_vals_combi);
        tb_pvals_combi.Properties.VariableNames = ls_st_param_key(cl_st_param_keys);
        tb_pvals_combi.Properties.RowNames = strcat('j=', string(1:size(mt_param_vals_combi,1)));
        clear mt_param_vals_combi;
        
        % Display
        disp(head(tb_pvals_combi, 10));
        disp(tail(tb_pvals_combi, 10));
        
        for it_pcombi_ctr = 1:it_test_combi_n
            
            tb_row_param_adj_cur = tb_pvals_combi(it_pcombi_ctr, :);
            
            disp(['xxxxx Shift to: xxxxx']);
            disp(tb_row_param_adj_cur)
            
            % Display Current Parameter been Varied
            for st_param_key = cl_st_param_keys
                param_map(st_param_key{1}) = tb_row_param_adj_cur{1, st_param_key};
            end
            
            % Simulate Model
            ar_simu_info = [it_pcombi_ctr, cell2mat(values(param_map, cl_st_param_keys))];
            cl_col_names = [{'it_pcombi_ctr'} cl_st_param_keys];
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
    
    % Table Output
    if (bl_display_simu_stats)
        disp('-------------------------');
        disp('xxxxxx tb_outcomes xxxxxx');
        disp(head(tb_outcomes, 10));
        disp(tail(tb_outcomes, 10));
    end
    
    %% Save Mat
    if (bl_mat_test)
        clear armt_map result_map
        if ~exist(st_mat_test_path,'dir'); mkdir(st_mat_test_path); end
        st_file_name = [st_mat_test_prefix st_mat_test_name_main st_mat_test_suffix];
        save(strcat(st_mat_test_path, st_file_name));
    end
    
end

end

%% Simulate Model given Parameter and Process Results
function [tb_outcomes_simu, variablenames] = simu_model_gen_stats(param_map, support_map, ar_simu_info, cl_col_names)

cl_st_param_keys = param_map('cl_st_param_keys');

%% Reset Parameters that are determined by other parameters
it_a_n = param_map('it_a_n');
it_z_n = param_map('it_z_n');
disp(['xxxxx it_a_n = ' num2str(it_a_n) ', it_z_n = ' num2str(it_z_n) ' xxxxx']);

%% Simulate Model
[armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map);
result_map = ff_az_vf_vecsv(param_map, support_map, armt_map, func_map);
result_map = ff_az_ds_vecsv(param_map, support_map, armt_map, func_map, result_map);

%% Store Results to Table

% Append Simulation Loop Statistics
tb_outcomes_meansdperc = result_map('tb_outcomes_meansdperc');
tb_simu_info = array2table(ar_simu_info + zeros([size(tb_outcomes_meansdperc,1), length(ar_simu_info)]));
tb_simu_info.Properties.VariableNames = cl_col_names;
variablenames = tb_outcomes_meansdperc.Properties.RowNames;
tb_simu_info.Properties.RowNames = (variablenames);

% Table
tb_outcomes_simu = [tb_simu_info tb_outcomes_meansdperc];
end
