%% Derive Savings Distribution (Loop)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_az_ds(varargin)
%% FF_AZ_DS finds the stationarz savings distribution
% Building on the Asset Dynamic Programming Problem
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html
% ff_az_vf_vecsv>, here we solve for the asset distribution. 
%
% This is the iterative looped solution to derive the asset distribution
% determined by the policy functions. Note that the asset distribution is a
% joint discrete random variable. We derive f(a,z), where f is the joint
% discrete random variables probability mass. Then we can derive f(a'(a,z)),
% f(c(a,z)) directly. The procedure here does not involve simulation.
% Simulation could also be used to derive these distributions, but given
% the discrete grid based solution algorithm, there is no need to
% introduce simulation and associated errors once we have fixed the shock
% process that generates randomness. 
%
% The function here accomplishes two tasks: (1) deriving the asset
% distribution as a discrete random variable over the states (2)
% calculating various statistics based on the discrete joint random
% variable's probability mass function for various outcomes of the model
%
% Distributions of Interest:
%
% * $p(a,z)$
% * $p(Y=y, z) = \sum_{a} \left( 1\left\{Y(a,z)=y\right\} \cdot p(a,z) \right)$
% * $p(Y=y) = \sum_{a,z} \left( 1\left\{Y(a,z)=y\right\} \cdot p(a,z) \right)$
%
% Statistics include:
%
% * $\mu_y = \sum_{y} p(Y=y) \cdot y$
% * $\sigma_y = \sqrt{ \sum_{y} p(Y=y) \cdot \left( y - \mu_y \right)^2}$
% * $p(y=0)$
% * $p(y=\max(y))$
% * percentiles: $min_{y} \left\{ P(Y \le y) - percentile \mid P(Y \le y) \ge percentile \right\}$
% * fraction of outcome held by up to percentiles: $E(Y<y)/E(Y)$
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param armt_map container container with states, choices and shocks
% grids that are inputs for grid based solution algorithm
%
% @param func_map container container with function handles for
% consumption cash-on-hand etc.
%
% @return result_map container contains policy function matrix, value
% function matrix, iteration results, and policy function, value function
% and iteration results tables.
%
% new keys included in result_map in addition to the output from
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html
% ff_az_vf_vecsv> are various distribution statistics for each model
% outcome, keys include *cl_mt_pol_a*, *cl_mt_pol_c*, *cl_mt_pol_coh*, etc,
% these include:
%
% * the first element of each of these cell array is y(a,z), the
% outcome/choice at the state space points
% * the second element of the cell is another container, which contains
% statistics computed for f(y) based on y(a,z) and f(a,z), f(y) is the
% probability mass function for outcome y given the stationary distribution
% f(a,z). The second element container also includes f(y) itself as well as
% f(y,z).
% * additionally, result_map also stores some of the statistics for
% different variables jointly together. (a) *tb_outcomes_meansdperc*: where
% each row is a different outcome of the model, and each table column
% stores a different statistics of interest. (b) *tb_outcomes_fracheld*:
% which measures the fraction of asset held by different people.
%
% @example
%
%    % Get Default Parameters
%    it_param_set = 6;
%    [param_map, support_map] = ffs_abz_set_default_param(it_param_set);
%    % Change Keys in param_map
%    param_map('it_a_n') = 500;
%    param_map('it_z_n') = 11;
%    param_map('fl_a_max') = 100;
%    param_map('fl_w') = 1.3;
%    % Change Keys support_map
%    support_map('bl_display') = false;
%    support_map('bl_post') = true;
%    support_map('bl_display_final') = false;
%    % Call Program with external parameters that override defaults
%    ff_az_ds(param_map, support_map);
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html ff_az_vf_vecsv>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html fft_disc_rand_var_stats>
%
% @seealso
%

%% Default
% # it_subset = 5 is basic invoke quick test
% # it_subset = 6 is main invoke
% # it_subset = 7 is profiling invoke
% # it_subset = 8 is matlab publish.

it_param_set = 8;
bl_input_override = true;
[param_map, support_map] = ffs_az_set_default_param(it_param_set);

% Note: param_map and support_map can be adjusted here or outside to override defaults
% param_map('it_a_n') = 750;
% param_map('it_z_n') = 15;

[armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
default_params = {param_map support_map armt_map func_map};

%% Parse Parameters 1

% if varargin only has param_map and support_map,
params_len = length(varargin);
[default_params{1:params_len}] = varargin{:};
param_map = [param_map; default_params{1}];
support_map = [support_map; default_params{2}];
if params_len >= 1 && params_len <= 2
    % If override param_map, re-generate armt and func if they are not
    % provided
    bl_input_override = true;
    [armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map, bl_input_override);
else
    % Override all
    armt_map = [armt_map; default_params{3}];
    func_map = [func_map; default_params{4}];
end

% append function name
st_func_name = 'ff_az_ds';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

%% Obtain Result Map By Solving Value/Pol Function
% use the optimized-vectorized ff_az_vf_vecsv.m function which works the
% fastest and produces identical results as ff_az_vf_vec and ff_az_vf. 

result_map = ff_az_vf_vecsv(param_map, support_map);

%% Parse Parameters 2

% result_map
% ar_st_pol_names is from section _Process Optimal Choices_ in the value
% function code.
params_group = values(result_map, {'cl_mt_pol_a', 'ar_st_pol_names'});
[cl_mt_pol_a, ar_st_pol_names] = params_group{:};
mt_pol_a = deal(cl_mt_pol_a{1});

% armt_map
params_group = values(armt_map, {'ar_a', 'mt_z_trans', 'ar_z'});
[ar_a, mt_z_trans, ar_z] = params_group{:};

% param_map
params_group = values(param_map, {'it_a_n', 'it_z_n'});
[it_a_n, it_z_n] = params_group{:};
params_group = values(param_map, {'it_maxiter_dist', 'fl_tol_dist'});
[it_maxiter_dist, fl_tol_dist] = params_group{:};

% support_map
params_group = values(support_map, {'bl_profile_dist', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time', 'bl_display_dist', 'it_display_every', 'bl_display_final_dist', 'bl_post'});
[bl_profile_dist, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display_dist, it_display_every, bl_display_final_dist, bl_post] = params_group{:};

%% Start Profiler and Timer

% Start Profile
if (bl_profile_dist)
    close all;
    profile off;
    profile on;
end

% Start Timer
if (bl_time)
    tic;
end

%% *f(a,z)*: Initialize Output Matrixes
% Initialize the distribution to be uniform

mt_dist_az_init = ones(length(ar_a),length(ar_z))/length(ar_a)/length(ar_z);
mt_dist_az_cur = mt_dist_az_init;
mt_dist_az_zeros = zeros(length(ar_a),length(ar_z));

%% *f(a,z)*: Initialize Convergence Conditions

bl_histiter_continue = true;
it_iter = 0;
ar_dist_diff_norm = zeros([it_maxiter_dist, 1]);
mt_dist_perc_change = zeros([it_maxiter_dist, it_z_n]);

%% *f(a,z)*: Derive Stationary Distribution
% Iterate over the discrete joint random variable variables (a,z)
while (bl_histiter_continue)
    
    it_iter = it_iter + 1;
    
    %% *f(a,z)*: Iterate over Probability mass for Discrete Random Variable
    % compared to
    % <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf.html
    % ff_az_vf>, we basically have the same set of loops. There, there were
    % four loops, here there are three loops. We eliminated the loop over
    % next period choices, because here we already know optimal choices
    
    % initialize empty
    mt_dist_az = mt_dist_az_zeros;
    
    % loop 1: over exogenous states
    for it_z_i = 1:length(ar_z)
        
        % loop 2: over endogenous states
        for it_a_j = 1:length(ar_a)
            
            % f(a'|a) = 1 for only one a'
            % in dynamic programming problem, had a loop over choices, now
            % already have optimal choices, do not need to loop
            fl_aprime = mt_pol_a(it_a_j, it_z_i);           
            it_a_prime_idx = find(ar_a == fl_aprime);
            
            % loop 3: loop over future shocks
            % E_{z'}(f(a',z'|a,z)*f(a,z))
            for it_zp_q = 1:length(ar_z)
                
                % current probablity at (a,z)
                fl_cur_za_prob = mt_dist_az_cur(it_a_j, it_z_i);
                
                % f(z'|z) transition
                fl_ztoz_trans =  mt_z_trans(it_z_i, it_zp_q);
                
                % f(a',z'|a,z)*f(a,z) 
                fl_zfromza = fl_cur_za_prob*fl_ztoz_trans;
                
                % cumulating
                mt_dist_az(it_a_prime_idx, it_zp_q) = mt_dist_az(it_a_prime_idx, it_zp_q) + fl_zfromza;
            end 
            
        end
        
    end
   
    
    %% *f(a,z)*: Check Tolerance and Continuation
    
    % Difference across iterations
    ar_dist_diff_norm(it_iter) = norm(mt_dist_az - mt_dist_az_cur);
    mt_dist_perc_change(it_iter, :) = sum((mt_dist_az ~= mt_dist_az))/(it_a_n);
    
    % Update
    mt_dist_az_cur = mt_dist_az;    
    
    % Print Iteration Results
    if (bl_display_dist && (rem(it_iter, it_display_every)==0))
        fprintf('Dist it_iter:%d, fl_dist_diff:%d\n', it_iter, ar_dist_diff_norm(it_iter));
        tb_hist_iter = array2table([sum(mt_dist_az_cur,1); std(mt_dist_az_cur,1); ...
                                    mt_dist_az_cur(1,:); mt_dist_az_cur(it_a_n,:)]);
        tb_hist_iter.Properties.VariableNames = strcat('z', string((1:size(mt_dist_az,2))));
        tb_hist_iter.Properties.RowNames = {'mdist','sddist', 'Ldist', 'Hdist'};
        disp('mdist = sum(mt_dist_az_cur,1) = sum_{a}(p(a)|z)')
        disp('sddist = std(mt_pol_a_cur,1) = std_{a}(p(a)|z)')
        disp('Ldist = mt_dist_az_cur(1,:) = p(min(a)|z)')
        disp('Hdist = mt_dist_az_cur(it_a_n,:) = p(max(a)|z)')
        disp(tb_hist_iter);
    end

    % Continuation Conditions:
    if (it_iter == (it_maxiter_dist + 1))
        bl_histiter_continue = false;
    elseif ((it_iter == it_maxiter_dist) || ...
            (ar_dist_diff_norm(it_iter) < fl_tol_dist))
        it_iter_last = it_iter;
        it_iter = it_maxiter_dist;
    end
    
end

%% End Time and Profiler

% End Timer
if (bl_time)
    toc;
end

% End Profile
if (bl_profile_dist)
    profile off
    profile viewer
    st_file_name = [st_profile_prefix st_profile_name_main st_profile_suffix];
    profsave(profile('info'), strcat(st_profile_path, st_file_name));
end

%% *f(y), f(c), f(a)*: Generate Key Distributional Statistics for Each outcome
% Having derived f(a,z) the probability mass function of the joint discrete
% random variables, we now obtain distributional statistics. Note that we
% know f(a,z), and we also know relevant policy functions a'(a,z), c(a,z),
% or other policy functions. We can simulate any choices that are a
% function of the random variables (a,z), using f(a,z)
%
% parameter structure provides a list of 
%
% # from result_map('ar_st_pol_names'), get list of outcome matrix on state
% space
% # simulate each outcome using f(a,z) for probability draws
% # compute key statistics: (1) mean (expectation=sum) (2) sd (3) min and
% max (4) iqr (5) fraction = 0 (6) percentiles including: 99.9, 99, 95,
% every 5 in between 5, 1, 0.01. 
%

% Loop over outcomes, see end of
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html
% ff_az_vf_vecsv> where these are created
for it_outcome_ctr=1:length(ar_st_pol_names)
    
    %% *f(y), f(c), f(a)*: Find p(outcome(states)), proability mass function for each outcome
    % Using from tools:
    % <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html
    % fft_disc_rand_var_mass2outcomes>, compute unique sorted outcomes for
    % y(a,z) and find:
    %
    % $$ p(y,z) = \sum_{a} \left(1\left\{Y(a,z)=y\right\} \cdot p(a,z) \right)$$
    %
    % $$ p(Y=y) = \sum_{a,z} \left( 1\left\{Y(a,z)=y\right\} \cdot p(a,z) \right)$$
    %
    % note: sum(mt_dist_az, 2) = result_map('cl_mt_pol_a'){2}, but not at
    % small simulation grids. These two might be different because pol_a is
    % based on a choices, mt_dist_az is based on a states
    %        
    % see end of
    % <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html
    % ff_az_vf_vecsv> outcomes in result_map are cells with two elements,
    % first element is y(a,z), second element will be f(y) and y, generated
    % here.
    %
    
    st_cur_output_key = ar_st_pol_names(it_outcome_ctr);
    cl_mt_choice_cur = result_map(st_cur_output_key);
    mt_choice_cur = cl_mt_choice_cur{1};
    
    % run function from tools: fft_disc_rand_var_mass2outcomes
    % <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html>
    bl_input_override = true;
    [tb_choice_drv_cur_byY, ar_choice_prob_byY, ar_choice_unique_sorted_byY, mt_choice_prob_byYZ] = ...
        fft_disc_rand_var_mass2outcomes(st_cur_output_key, mt_choice_cur, mt_dist_az, bl_input_override);
    
    %% *f(y), f(c), f(a)*: Compute Statistics for outcomes
    % Using from tools:
    % <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html
    % fft_disc_rand_var_stats>, compute these outcomes:
    %
    % * $\mu_Y = E(Y) = \sum_{y} p(Y=y) \cdot y $
    % * $\sigma_Y = \sqrt{ \sum_{y} p(Y=y) \cdot \left( y - \mu_y \right)^2}$
    % * $p(y=0)$
    % * $p(y=\max(y))$
    % * percentiles: $min_{y} \left\{ P(Y \le y) - percentile \mid P(Y \le y) \ge percentile \right\}$
    % * fraction of outcome held by up to percentiles: $E(Y<y)/E(Y)$
    %
    
    % run function fft_disc_rand_var_stats.m from tools:
    % <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html>
    [ds_stats_map] = fft_disc_rand_var_stats(st_cur_output_key, ar_choice_unique_sorted_byY', ar_choice_prob_byY');
    
    % prcess results
    % retrieve scalar statistics: 
    fl_choice_mean = ds_stats_map('fl_choice_mean');
    fl_choice_sd = ds_stats_map('fl_choice_sd');
    fl_choice_coefofvar = ds_stats_map('fl_choice_coefofvar');
    fl_choice_min = ds_stats_map('fl_choice_min');
    fl_choice_max = ds_stats_map('fl_choice_max');
    fl_choice_prob_zero = ds_stats_map('fl_choice_prob_zero');
    fl_choice_prob_min = ds_stats_map('fl_choice_prob_min');
    fl_choice_prob_max = ds_stats_map('fl_choice_prob_max');
    % retrieve distributional array stats
    tb_prob_drv = ds_stats_map('tb_prob_drv');
    ar_choice_percentiles = tb_prob_drv{:,2};
    ar_choice_perc_fracheld = tb_prob_drv{:,3};   

    % Display
    if (bl_display_final_dist)
        disp(['tb_prob_drv, Percentiles of Y, and Share of Y Held by Households up to this Percentile: ', st_cur_output_key])
        disp(tb_prob_drv);
    end
    
    %% *f(y), f(c), f(a)*: Store Statistics Specific to Each Outcome
    % see intro section
    
    % Append prob mass functions to ds_stats_map
    ds_stats_map('mt_choice_prob_byYZ') = mt_choice_prob_byYZ;
    ds_stats_map('tb_choice_prob_byY') = tb_choice_drv_cur_byY;
    % ds_stats_map is second element of cell for the key for the variable
    % in result_map
    cl_mt_choice_cur{2} = ds_stats_map;    
    result_map(st_cur_output_key) = cl_mt_choice_cur;
    
    % key stats
    ar_keystats = [fl_choice_mean fl_choice_sd fl_choice_coefofvar fl_choice_min fl_choice_max ...
        fl_choice_prob_zero fl_choice_prob_min fl_choice_prob_max ar_choice_percentiles'];
    cl_outcome_names(it_outcome_ctr) = st_cur_output_key;
    if (it_outcome_ctr == 1)
        mt_outcomes_meansdperc = ar_keystats;
        mt_outcomes_fracheld = ar_choice_perc_fracheld';
    else
        mt_outcomes_meansdperc = [mt_outcomes_meansdperc; ar_keystats];
        mt_outcomes_fracheld = [mt_outcomes_fracheld; ar_choice_perc_fracheld'];
    end
    
end

% *f(y), f(c), f(a)*: Store Statistics Shared Table All Outcomes
% Process mean and and percentiles
tb_outcomes_meansdperc = array2table(mt_outcomes_meansdperc);
ar_fl_percentiles = tb_prob_drv{:,1};
cl_col_names = ['mean', 'sd', 'coefofvar', 'min', 'max', ...
                'pYis0', 'pYisMINY', 'pYisMAXY', strcat('p', string(ar_fl_percentiles'))];
tb_outcomes_meansdperc.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
tb_outcomes_meansdperc.Properties.RowNames = matlab.lang.makeValidName(cl_outcome_names);

% Process Aset Held by up to percentiles
tb_outcomes_fracheld = array2table(mt_outcomes_fracheld);
ar_fl_percentiles = tb_prob_drv{:,1};
cl_col_names = [strcat('fracByP', string(ar_fl_percentiles'))];
tb_outcomes_fracheld.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
tb_outcomes_fracheld.Properties.RowNames = matlab.lang.makeValidName(cl_outcome_names);

% Add to result_map
result_map('tb_outcomes_meansdperc') = tb_outcomes_meansdperc;
result_map('mt_outcomes_fracheld') = mt_outcomes_fracheld;

% Display
if (bl_display_final_dist)
    
    disp('tb_outcomes_meansdperc: mean, sd, percentiles')
    disp(tb_outcomes_meansdperc);
    
    disp('tb_outcomes_fracheld: fraction of asset/income/etc held by hh up to this percentile')
    disp(tb_outcomes_fracheld);
        
end


%% Process Optimal Choices

% 
% if (bl_post)
%     bl_input_override = true;
%     result_map('ar_dist_diff_norm') = ar_dist_diff_norm(1:it_iter_last);
%     result_map('mt_dist_perc_change') = mt_dist_perc_change(1:it_iter_last, :);
%     result_map = ff_ds_vf_post(param_map, support_map, armt_map, func_map, result_map, bl_input_override);
% end
% 
end
