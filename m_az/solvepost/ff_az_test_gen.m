%% Simulate Savings Model Along Various Parameter Arrays
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [tb_outcomes, support_map, param_desc_map] = ff_az_test_gen(varargin)
%% FF_AZ_TEST_GEN post solution simulation
% Simulate the model along various dimensions
%
% @param st_simu_type string cross vs gridd simulation. cross: with
% (x,y), vary X fixing y, vary Y fixing x. grid: all (x,y) \in (X,Y)
%
% # _st_simu_type_ = 'c' for cross simulate, if 'c', then each array value
% of param_tstar_map container is simulated over one by one. So if there
% are two arrays associated with two keys in param_tstar_map with length N1
% and N2, the total number of simulations equals to N1 + N2.
% # _st_simu_type_ = 'g' for grid simulate, if 'g', then all possible
% combinations of the arrays in param_tstar_map are used to create
% combinations of parameter values. So if there are two arrays with lengths
% N1 and N2, there will be N1*N2 number of simulations
% # _st_simu_type_ = 'r' for random simulate, if 'r', then should specify
% param_map('it_st_simu_type_g_seed') value. Only the minimum and the
% maximum values for each array in param_tstar_map matters. Based on these
% minimum and maximum, and also what is in
% param_map('it_st_simu_type_g_simun'). Random parameter values will be
% drawn.
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
st_simu_type = 'r'; % 'c' for cross or 'g' for grid or 'r' rand

%% Default Parameters

% Base parameters call
it_param_set = 9;
[param_map, support_map] = ffs_az_set_default_param(it_param_set);

% allows verticle cleaner display of stats. Which is useful.
support_map('bl_display_final_dist_detail') = true;

% Timer
support_map('bl_timer') = true;

% Mat Storage
support_map('bl_mat_test') = true;
st_matimg_path_root = support_map('st_matimg_path_root');
% test_borinf for default: cl_st_param_keys = [3, 8, 9]
support_map('st_mat_test_path') = [st_matimg_path_root '/test/ff_az_ds_vecsv/mat/'];
support_map('st_mat_test_prefix') = [''];

% Generate or Load
% if bl_replacefile is true, that means even if file already exists,
% generaet new. bl_replacefile is false, then if false does not exist
% generate, if file exists does not generate.
support_map('bl_replacefile') = true;
support_map('bl_display_simu_stats') = true;

% param map modification
param_map('it_st_simu_type_g_seed') = 123;
param_map('it_st_simu_type_g_simun') = 20;

%% Array Parameters
% keys here for az, abz models

ls_st_param_key = {'fl_crra', 'fl_beta', ...
                    'fl_w', 'fl_r_save', ...
                    'fl_z_rho', 'fl_z_sig', ...
                    'fl_a_max', 'it_z_n', 'it_a_n', ...
                    ...
                    'fl_z_r_borr_poiss_mean', 'fl_z_r_borr_max', ...
                    'fl_b_bd', 'fl_c_min', 'fl_z_r_borr_n',...
                    'fl_z_wage_rho', 'fl_z_wage_sig', ...
                    ...
                    'fl_alpha', 'fl_delta', ...
                    ...
                    'fl_r_borr'};

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

param_tstar_map(ls_st_param_key{10}) = linspace(5, 20, it_simu_vec_len);
param_tstar_map(ls_st_param_key{11}) = linspace(0.095, 0.150, it_simu_vec_len);
param_tstar_map(ls_st_param_key{12}) = linspace(-20, -5, it_simu_vec_len);
param_tstar_map(ls_st_param_key{13}) = linspace(0.03, 0.001, it_simu_vec_len);
param_tstar_map(ls_st_param_key{14}) = 5:4:25;
param_tstar_map(ls_st_param_key{15}) = linspace(0, 0.99, it_simu_vec_len);
param_tstar_map(ls_st_param_key{16}) = linspace(0.01, 0.5, it_simu_vec_len);

param_tstar_map(ls_st_param_key{17}) = linspace(0.30, 0.50, it_simu_vec_len);
param_tstar_map(ls_st_param_key{18}) = linspace(0.02, 0.14, it_simu_vec_len);

param_tstar_map(ls_st_param_key{19}) = linspace(0, 0.25, it_simu_vec_len);

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

param_desc_map(ls_st_param_key{10}) = {'pois-max(shift br r mean)'};
param_desc_map(ls_st_param_key{11}) = {'max borrow r'};
param_desc_map(ls_st_param_key{12}) = {'borrow bound'};
param_desc_map(ls_st_param_key{13}) = {'minimum consumption'};
param_desc_map(ls_st_param_key{14}) = {'borrow shock n'};
param_desc_map(ls_st_param_key{15}) = {'Inc/Prod Shock Persistence'};
param_desc_map(ls_st_param_key{16}) = {'Inc/Prod Shock SD'};

param_desc_map(ls_st_param_key{17}) = {'Prod Func Elasticity'};
param_desc_map(ls_st_param_key{18}) = {'Prod Func Depreciation'};

param_desc_map(ls_st_param_key{19}) = {'Single Borrow Rate'};

%% Default

default_params = {st_simu_type it_size_type cl_st_param_keys ...
                    param_map support_map param_tstar_map param_desc_map};

%% Parse Inputs
% override default set above if any parameters are updated

params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};

st_simu_type = default_params{1};
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

params_group = values(param_map, {'st_model'});
[st_model] = params_group{:};

params_group = values(support_map, {'bl_replacefile', 'bl_display_simu_stats'});
[bl_replacefile, bl_display_simu_stats] = params_group{:};

params_group = values(param_map, {'it_st_simu_type_g_seed', 'it_st_simu_type_g_simun'});
[it_st_simu_type_g_seed, it_st_simu_type_g_simun] = params_group{:};

%% Store Mat Strings
%
% # _it_total_length_: this is the total length of all selected arrays
% # _ar_param_keys_idx_: these are the linear index of the selected keys
% among all possible keys
%

it_total_length = sum(cell2mat(cellfun(@(m) length(param_tstar_map(m)), ...
                                cl_st_param_keys, 'UniformOutput', false)));

support_map('st_mat_test_name_main') = ['r'];
support_map('st_mat_test_suffix') = ['_g' strrep(num2str(ar_param_keys_idx), '  ', '') ...
                                     '_c' st_simu_type ...
                                     '_t' num2str(it_size_type) 'l' num2str(it_total_length)];

%% Set Solve Sizes

if (it_size_type == 1)

    % Basic Test Run
    if (ismember(st_model, ["az"]))

        param_map('it_z_n') = 11;

        param_map('it_a_n') = 100;

    elseif (ismember(st_model, ["abz", "ipwkbzr"]))

        if (ismember(st_model, ["abz"]))

            param_map('fl_z_r_borr_n') = 5;
            param_map('it_z_wage_n') = 5;
            param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

            param_map('it_a_n') = 100;

        elseif (ismember(st_model, ["ipwkbzr"]))

            param_map('fl_z_r_borr_n') = 5;
            param_map('it_z_wage_n') = 5;
            param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');

            param_map('fl_coh_interp_grid_gap') = 0.3;
            param_map('it_c_interp_grid_gap') = 10^-4;
            param_map('it_w_perc_n') = 25;
            param_map('it_ak_perc_n') = param_map('it_w_perc_n');
            param_map('fl_w_interp_grid_gap') = 0.3;

        end

    elseif (ismember(st_model, ["akz_wkz_iwkz", "ipwkz"]))

        param_map('it_z_n') = 11;
        param_map('fl_coh_interp_grid_gap') = 0.2;
        param_map('it_c_interp_grid_gap') = 10^-4;

        if (ismember(st_model, ["akz_wkz_iwkz"]))

            param_map('it_w_n') = 25;
            param_map('it_ak_n') = param_map('it_w_n');

        elseif (ismember(st_model, ["ipwkz"]))

            param_map('it_w_perc_n') = 25;
            param_map('it_ak_perc_n') = param_map('it_w_perc_n');
            param_map('fl_w_interp_grid_gap') = 0.2;

        end

    end

elseif (it_size_type == 2)

    % Full Run

elseif (it_size_type == 3)

    % Denser Run
    if (ismember(st_model, ["az"]))

        param_map('it_z_n') = 27;

        param_map('it_a_n') = 2250;

    elseif (ismember(st_model, ["abz", "ipwkbzr"]))


        if (ismember(st_model, ["abz"]))

            % keep the interest rate structure the same as default
            param_map('fl_z_r_borr_n') = 11;
            param_map('it_z_wage_n') = 11;
            param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
            
            param_map('it_a_n') = 1250;

        elseif (ismember(st_model, ["ipwkbzr"]))

            % keep the interest rate structure the same as default
            param_map('fl_z_r_borr_n') = 7;
            param_map('it_z_wage_n') = 11;
            param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
            
            param_map('fl_coh_interp_grid_gap') = 0.07;
            param_map('it_c_interp_grid_gap') = 10^-4;
            param_map('it_w_perc_n') = 85;
            param_map('it_ak_perc_n') = param_map('it_w_perc_n');
            param_map('fl_w_interp_grid_gap') = 0.07;

        end

    elseif (ismember(st_model, ["akz_wkz_iwkz", "ipwkz"]))

        param_map('it_z_n') = 21;
        param_map('fl_coh_interp_grid_gap') = 0.025;
        param_map('it_c_interp_grid_gap') = 10^-4;

        if (ismember(st_model, ["akz_wkz_iwkz"]))

            param_map('it_w_n') = 150;
            param_map('it_ak_n') = param_map('it_w_n');

        elseif (ismember(st_model, ["ipwkz"]))

            param_map('it_w_perc_n') = 150;
            param_map('it_ak_perc_n') = param_map('it_w_perc_n');
            param_map('fl_w_interp_grid_gap') = 0.025;

        end

    end
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

    if (strcmp(st_simu_type, 'c'))

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
            for it_param_all_ctr = 1:length(cl_st_param_keys)
                disp([ cl_st_param_keys{it_param_all_ctr} ':' num2str(param_map(cl_st_param_keys{it_param_all_ctr}))]);
            end
            disp('xxxxxxxxxxxxxxxxxxxxxxxxxxx');

            % Simulate Model over Parameter Array Values
            for it_cur_param = 1:1:length(ar_param_values)

                % Adjust Value for Current Parameter been Varied
                fl_param_val = ar_param_values(it_cur_param);
                param_map(st_param_key) = fl_param_val;
                param_map = map_correction(param_map);
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

    if (ismember(st_simu_type, ["g", "r"]))

        % Get Arrays to be Meshed in Cells
        cl_ar_param_subset_values = values(param_tstar_map, cl_st_param_keys);

        % Generate all possible combinations of parameters subsets
        if (strcmp(st_simu_type, 'g'))
            cl_mt_all = cl_ar_param_subset_values;
            [cl_mt_all{:}] = ndgrid(cl_ar_param_subset_values{:});
            mt_param_vals_combi = cell2mat(cellfun(@(m) m(:), cl_mt_all, 'uni', 0));
        elseif (strcmp(st_simu_type, 'r'))
            % random draw within max and min N count
            rng(it_st_simu_type_g_seed);
            mt_param_vals_combi = cell2mat(cellfun(@(m) ...
                                       rand([it_st_simu_type_g_simun,1]).*(max(param_tstar_map(m)) - min(param_tstar_map(m))) ...
                                       + min(param_tstar_map(m)), ...
                                       cl_st_param_keys, 'UniformOutput', false));
        end

        % Sizes
        it_test_combi_n = size(mt_param_vals_combi,1);

        % Show Combinations of Parameters Simulating over, convert from Matrix to Table for Clarity
        tb_pvals_combi = array2table(mt_param_vals_combi);
        tb_pvals_combi.Properties.VariableNames = cl_st_param_keys;
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
            param_map = map_correction(param_map);

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
st_model = param_map('st_model');
if (ismember(st_model, ["az","abz"]))
    it_a_n = param_map('it_a_n');
    it_z_n = param_map('it_z_n');
    disp(['xxxxx st_model = ' st_model ', it_a_n = ' num2str(it_a_n) ', it_z_n = ' num2str(it_z_n) ' xxxxx']);
elseif (ismember(st_model, ["akz_wkz_iwkz"]))
    it_w_n = param_map('it_w_n');
    it_ak_n = param_map('it_ak_n');
    it_z_n = param_map('it_z_n');
    disp(['xxxxx st_model = ' st_model ...
          ', it_w_n = ' num2str(it_w_n) ', it_ak_n = ' num2str(it_ak_n) ...
          ', it_z_n = ' num2str(it_z_n) ' xxxxx']);
elseif (ismember(st_model, ["ipwkz", "ipwkbzr"]))
    it_w_perc_n = param_map('it_w_perc_n');
    it_ak_perc_n = param_map('it_ak_perc_n');
    it_z_n = param_map('it_z_n');
    disp(['xxxxx st_model = ' st_model ...
          ', it_w_perc_n = ' num2str(it_w_perc_n) ', it_ak_perc_n = ' num2str(it_ak_perc_n) ...
          ', it_z_n = ' num2str(it_z_n) ' xxxxx']);
end


%% Simulate Model

if (ismember(st_model, ["az", "abz"]))

    if (ismember(st_model, ["az"]))
        [armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map);
        result_map = ff_az_vf_vecsv(param_map, support_map, armt_map, func_map);
    elseif (ismember(st_model, ["abz"]))
        [armt_map, func_map] = ffs_abz_get_funcgrid(param_map, support_map);
        result_map = ff_abz_vf_vecsv(param_map, support_map, armt_map, func_map);
    end

    result_map = ff_az_ds_vec(param_map, support_map, armt_map, func_map, result_map);

elseif (ismember(st_model, ["akz_wkz_iwkz"]))

    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map);
    result_map = ff_iwkz_vf_vecsv(param_map, support_map, armt_map, func_map);
    result_map = ff_iwkz_ds_vec(param_map, support_map, armt_map, func_map, result_map);

elseif (ismember(st_model, ["ipwkz"]))

    [armt_map, func_map] = ffs_ipwkz_get_funcgrid(param_map, support_map);
    result_map = ff_ipwkz_vf_vecsv(param_map, support_map, armt_map, func_map);
    result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map);

elseif (ismember(st_model, ["ipwkbzr"]))

    [armt_map, func_map] = ffs_ipwkbzr_get_funcgrid(param_map, support_map);
    result_map = ff_ipwkbzr_vf_vecsv(param_map, support_map, armt_map, func_map);
    result_map = ff_iwkz_ds_vecsv(param_map, support_map, armt_map, func_map, result_map);

end

%% Store Results to Table

% Append Simulation Loop Statistics
tb_outcomes_meansdperc = result_map('tb_outcomes');
tb_simu_info = array2table(ar_simu_info + zeros([size(tb_outcomes_meansdperc,1), length(ar_simu_info)]));
tb_simu_info.Properties.VariableNames = cl_col_names;
variablenames = tb_outcomes_meansdperc.Properties.RowNames;
tb_simu_info.Properties.RowNames = (variablenames);

% Table
tb_outcomes_simu = [tb_simu_info tb_outcomes_meansdperc];
end

%% Map Corrections
function param_map = map_correction(param_map)

    params_group = values(param_map, {'st_model'});
    [st_model] = params_group{:};
    
    if (isKey(param_map, 'fl_z_r_borr_n'))
        param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
    end

    if (ismember(st_model, ["ipwkbz", "ipwkbzr"]))

        if (isKey(param_map, 'fl_b_bd'))
            param_map('fl_w_min') = param_map('fl_b_bd');
            param_map('fl_k_max') = (param_map('fl_w_max') - param_map('fl_b_bd'));
        end
        if (isKey(param_map, 'fl_w_max'))
            param_map('fl_k_max') = (param_map('fl_w_max') - param_map('fl_b_bd'));
        end
    end
    
    if (ismember(st_model, ["abz", "ipwkbzr"]))    
        % Set Borrowing and Savings Interest Rates the same
        
        if (isKey(param_map, 'fl_r_borr'))
            % This key does not have any effect on abz or ipwkbzr, the
            % intent is to generate a single borrowing interest rate.
            param_map('fl_z_r_borr_max') = param_map('fl_r_borr');
            param_map('fl_z_r_borr_min') = param_map('fl_r_borr');
            param_map('fl_z_r_borr_n') = 1;            
            param_map('it_z_n') = param_map('it_z_wage_n') * param_map('fl_z_r_borr_n');
        end
    end


end
