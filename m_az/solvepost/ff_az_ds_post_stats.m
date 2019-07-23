%% Generate Statistics from Prob Mass Function over States
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_az_ds_post_stats(varargin)
%% FF_AZ_DS_POST_STATS post ff_az_ds statistics generation 
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
% Uses fake binomial data when file is invoke with defaults.
%
% @param param_map container parameter container
%
% @param support_map container support container
%
% @param result_map container contains policy function matrix, value
% function matrix, iteration results
%
% @param mt_dist_az matrix N by M where N are asset states and M are shock
% states, the f(a,z) probability mass function derived earlier in ff_az_ds
% or ff_az_ds_vec
%
% @return result_map container with statistics added to result_map
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
%    bl_input_override = true;
%    result_map = ff_az_ds_post_stats(support_map, result_map, mt_dist_az, bl_input_override);
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html fft_disc_rand_var_stats>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html fft_disc_rand_var_mass2outcomes>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2covcor.html fft_disc_rand_var_mass2covcor>
%

%% Default
% use binomial as test case, z maps to binomial win prob, remember binom
% approximates normal.

params_len = length(varargin);
bl_input_override = 0;
if (params_len == 4)
    bl_input_override = varargin{4};
end

if (bl_input_override)
    % if invoked from outside overrid fully
    [support_map, result_map, mt_dist_az, ~] = varargin{:};
    bl_display_final_dist_detail_local = false;
else
    clear all;
    close all;
    
    it_states = 6;
    it_shocks = 5;
    fl_binom_n = it_states-1;
    ar_binom_p = (1:(it_shocks))./(it_shocks+2);
    ar_binom_x = (0:1:(it_states-1)) -3;
    
    % f(z)
    ar_binom_p_prob = binopdf(0:(it_shocks-1), it_shocks-1, 0.5);
    % f(a,z), mass for a, z
    mt_dist_az = zeros([it_states, it_shocks]);
    for it_z=1:it_shocks
        % f(a|z)
        f_a_condi_z = binopdf(ar_binom_x - min(ar_binom_x), fl_binom_n, ar_binom_p(it_z));
        % f(z)
        f_z = ar_binom_p_prob(it_z);
        % f(a,z)=f(a|z)*f(z)
        mt_dist_az(:, it_z) = f_a_condi_z*f_z;
    end
    
    % y(a,z), some non-smooth structure
    rng(123);    
    mt_pol_a = ar_binom_x' - 0.01*ar_binom_x'.^2  + ar_binom_p - 0.5*ar_binom_p.^2 + rand([it_states, it_shocks]);
    mt_pol_a = round(mt_pol_a*3);

    rng(456);
    mt_pol_c = 10 -(mt_pol_a) + 15*(rand([it_states, it_shocks])-0.5);
        
    % Generate result_map
    result_map = containers.Map('KeyType','char', 'ValueType','any');
    result_map('cl_mt_pol_a') = {mt_pol_a, zeros(1)};
    result_map('cl_mt_pol_c') = {mt_pol_c, zeros(1)};
    result_map('ar_st_pol_names') = ["cl_mt_pol_a", "cl_mt_pol_c"];
    
    % support_map
    support_map = containers.Map('KeyType','char', 'ValueType','any');    
    support_map('bl_display_final_dist') = true;
    support_map('bl_display_final_dist_detail') = true;
    bl_display_final_dist_detail_local = true;
end

%% Parse

% support_map
params_group = values(support_map, {'bl_display_final_dist', 'bl_display_final_dist_detail'});
[bl_display_final_dist, bl_display_final_dist_detail] = params_group{:};
if (bl_display_final_dist_detail)
    bl_display_drvstats = true;
else
    bl_display_drvstats = false;
end

% result_map
params_group = values(result_map, {'ar_st_pol_names'});
[ar_st_pol_names] = params_group{:};

%% *f(y), f(c), f(a)*: Generate Key Distributional Statistics for Each outcome
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
    % $$ p(y,a) = \sum_{z} \left(1\left\{Y(a,z)=y\right\} \cdot p(a,z) \right)$$
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
    [ar_choice_prob_byY, ar_choice_unique_sorted_byY, mt_choice_prob_byYZ, mt_choice_prob_byYA] = ...
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
    [ds_stats_map] = fft_disc_rand_var_stats(st_cur_output_key, ar_choice_unique_sorted_byY', ar_choice_prob_byY', bl_display_drvstats);

    % prcess results
    % retrieve scalar statistics:
    fl_choice_mean = ds_stats_map('fl_choice_mean');
    fl_choice_sd = ds_stats_map('fl_choice_sd');
    fl_choice_coefofvar = ds_stats_map('fl_choice_coefofvar');
    fl_choice_min = ds_stats_map('fl_choice_min');
    fl_choice_max = ds_stats_map('fl_choice_max');    
    fl_choice_prob_zero = ds_stats_map('fl_choice_prob_zero');
    fl_choice_prob_below_zero = ds_stats_map('fl_choice_prob_below_zero');
    fl_choice_prob_above_zero = ds_stats_map('fl_choice_prob_above_zero');        
    fl_choice_prob_min = ds_stats_map('fl_choice_prob_min');
    fl_choice_prob_max = ds_stats_map('fl_choice_prob_max');
    % retrieve distributional array stats
    ar_choice_percentiles = ds_stats_map('ar_choice_percentiles');
    ar_choice_perc_fracheld = ds_stats_map('ar_choice_perc_fracheld');

    %% *f(y), f(c), f(a)*: Store Statistics Specific to Each Outcome
    % see intro section

    % Append prob mass functions to ds_stats_map
    ds_stats_map('mt_choice_prob_byYZ') = mt_choice_prob_byYZ;
    ds_stats_map('mt_choice_prob_byYA') = mt_choice_prob_byYA;
    ds_stats_map('ar_choice_unique_sorted_byY') = ar_choice_unique_sorted_byY;
    ds_stats_map('ar_choice_prob_byY') = ar_choice_prob_byY;
    % ds_stats_map is second element of cell for the key for the variable
    % in result_map
    cl_mt_choice_cur{2} = ds_stats_map;
    result_map(st_cur_output_key) = cl_mt_choice_cur;

    % key stats
    ar_keystats = [fl_choice_mean fl_choice_sd fl_choice_coefofvar fl_choice_min fl_choice_max ...
        fl_choice_prob_zero fl_choice_prob_below_zero fl_choice_prob_above_zero ...
        fl_choice_prob_min fl_choice_prob_max ar_choice_percentiles];
    cl_outcome_names(it_outcome_ctr) = st_cur_output_key;
    if (it_outcome_ctr == 1)
        mt_outcomes_meansdperc = ar_keystats;
        mt_outcomes_fracheld = ar_choice_perc_fracheld;
    else
        mt_outcomes_meansdperc = [mt_outcomes_meansdperc; ar_keystats];
        mt_outcomes_fracheld = [mt_outcomes_fracheld; ar_choice_perc_fracheld];
    end

end

if (bl_display_final_dist || bl_display_final_dist_detail)
    
    disp('xxx PERCENTILES AND STATS xxx')        
    tb_outcomes_meansdperc = array2table(mt_outcomes_meansdperc);
    ar_fl_percentiles = ds_stats_map('ar_fl_percentiles');
    cl_col_names = ['mean', 'sd', 'coefofvar', 'min', 'max', ...
                    'pYis0', 'pYls0', 'pYgr0', 'pYisMINY', 'pYisMAXY', ...
                    strcat('p', string(ar_fl_percentiles))];
    tb_outcomes_meansdperc.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_outcomes_meansdperc.Properties.RowNames = matlab.lang.makeValidName(cl_outcome_names);
    disp('tb_outcomes_meansdperc: mean, sd, percentiles')
    
    if (bl_display_final_dist_detail_local)    
        disp(rows2vars(tb_outcomes_meansdperc));
    end

    % Process Aset Held by up to percentiles    
    disp('xxx All Variables Fraction of Y Held up to Percentile xxx')
    tb_outcomes_fracheld = array2table(mt_outcomes_fracheld);
    cl_col_names = [strcat('fracByP', string(ar_fl_percentiles))];
    tb_outcomes_fracheld.Properties.VariableNames = matlab.lang.makeValidName(cl_col_names);
    tb_outcomes_fracheld.Properties.RowNames = matlab.lang.makeValidName(cl_outcome_names);
    disp('tb_outcomes_fracheld: fraction of asset/income/etc held by hh up to this percentile')
    
    if (bl_display_final_dist_detail_local)    
        disp(rows2vars(tb_outcomes_fracheld));
    end
    
end

%% Covariance and Correlation
% Having computed elsewhere E(X), E(Y), and SD(X), SD(Y), and given X(a,z)
% and Y(a,z), which are the optimal choices along the endogenous state
% space grid a, and the exogenous state space grid z, and given also
% f(a,z), the probability mass function over (a,z), we compute covariance
% and correlation between outcomes X and Y. 
%
% * Covariance
%
% $$\mathrm{Cov}\left(x,y\right) = \sum_{a} \sum_{z} f(a,z) \cdot \left( x(a,z) - \mu_x \right) \cdot \left( y(a,z) - \mu_y \right)$$
%
% * Correlation
%
% $$\rho_{x,y} = \frac{\mathrm{Cov}\left(x,y\right)}{\sigma_x \cdot \sigma_y}$$

for it_outcome_x_ctr=1:length(ar_st_pol_names)
        
    st_cur_output_x_key = ar_st_pol_names(it_outcome_x_ctr);
    
    cl_mt_choice_cur = result_map(st_cur_output_x_key);
    ds_stats_map = cl_mt_choice_cur{2};

    cl_mt_choice_cur = result_map(st_cur_output_x_key);
    mt_choice_x_bystates = cl_mt_choice_cur{1};    
    fl_choice_x_mean = ds_stats_map('fl_choice_mean');
    fl_choice_x_sd = ds_stats_map('fl_choice_sd');    
    
    ar_covvar = zeros([1,length(ar_st_pol_names)*2]);
    ar_st_covvar = strings([1,length(ar_st_pol_names)*2]);
    for it_outcome_y_ctr=1:length(ar_st_pol_names)

        st_cur_output_y_key = ar_st_pol_names(it_outcome_y_ctr);

        cl_mt_choice_cur = result_map(st_cur_output_y_key);
        ds_stats_map = cl_mt_choice_cur{2};

        cl_mt_choice_cur = result_map(st_cur_output_y_key);
        mt_choice_y_bystates = cl_mt_choice_cur{1};        
        fl_choice_y_mean = ds_stats_map('fl_choice_mean');
        fl_choice_y_sd = ds_stats_map('fl_choice_sd');
        
        covvar_input_map = containers.Map('KeyType','char', 'ValueType','any');
        covvar_input_map('mt_choice_x_bystates') = mt_choice_x_bystates;
        covvar_input_map('mt_choice_y_bystates') = mt_choice_y_bystates;
        covvar_input_map('mt_dist_bystates') = mt_dist_az;
        covvar_input_map('fl_choice_x_mean') = fl_choice_x_mean;
        covvar_input_map('fl_choice_x_sd') = fl_choice_x_sd;
        covvar_input_map('fl_choice_y_mean') = fl_choice_y_mean;
        covvar_input_map('fl_choice_y_sd') = fl_choice_y_sd;
        
        [fl_cov_xy, fl_cor_xy] = fft_disc_rand_var_mass2covcor(covvar_input_map);
        
        % only include the y name, x name is from the row
        st_x_y_cov = strjoin(["fl_cov_" st_cur_output_y_key], '');
        st_x_y_cor = strjoin(["fl_cor_" st_cur_output_y_key], '');
        ds_stats_map(st_x_y_cov) = fl_cov_xy;
        ds_stats_map(st_x_y_cor) = fl_cor_xy;
        
        ar_covvar(it_outcome_y_ctr*2-1) = fl_cov_xy;
        ar_covvar(it_outcome_y_ctr*2) = fl_cor_xy;
        ar_st_covvar(it_outcome_y_ctr*2-1) = string(st_x_y_cov);
        ar_st_covvar(it_outcome_y_ctr*2) = string(st_x_y_cor);
        
        cl_mt_choice_cur{2} = ds_stats_map;
        result_map(st_cur_output_y_key) = cl_mt_choice_cur;        
    end    
    
    if (it_outcome_x_ctr == 1)
        mt_outcomes_covvar = ar_covvar;
    else
        mt_outcomes_covvar = [mt_outcomes_covvar; ar_covvar];
    end
    
end

if (bl_display_final_dist || bl_display_final_dist_detail)
    
    disp('xxx Variance and Covariance xxx')    
    tb_outcomes_covvar = array2table(mt_outcomes_covvar);
    tb_outcomes_covvar.Properties.VariableNames = matlab.lang.makeValidName(ar_st_covvar);
    tb_outcomes_covvar.Properties.RowNames = matlab.lang.makeValidName(cl_outcome_names);
    disp('tb_outcomes_covvar: variance correlation')
    
    if (bl_display_final_dist_detail_local)    
        disp(rows2vars(tb_outcomes_covvar));
    end
    
end

%% *f(y), f(c), f(a)*: Store Statistics Shared Table All Outcomes

% Add to result_map
mt_outcomes = [mt_outcomes_meansdperc, mt_outcomes_covvar, mt_outcomes_fracheld];
result_map('mt_outcomes') = mt_outcomes;

if (bl_display_final_dist || bl_display_final_dist_detail)
    
    tb_outcomes = [tb_outcomes_meansdperc, tb_outcomes_covvar, tb_outcomes_fracheld];
    result_map('tb_outcomes') = tb_outcomes;        
    
    disp(rows2vars(tb_outcomes));
    
end
end
