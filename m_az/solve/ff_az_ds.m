%% Derive Asset and Choices/Outcomes Distribution (Loop)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_az_ds(varargin)
%% FF_AZ_DS finds the stationary asset distributions
% Building on the Asset Dynamic Programming Problem
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html
% ff_az_vf_vecsv>, here we solve for the asset distribution. This version
% of the program uses loops.
%
% This finds the asset distribution induced by the policy functions. Note
% that the asset distribution is a joint discrete random variable. We
% derive f(a,z), where f is the joint discrete random variables probability
% mass. Then we can derive f(a'(a,z)), f(c(a,z)) directly. The procedure
% here does not involve simulation. Simulation could also be used to derive
% these distributions, but given the discrete grid based solution
% algorithm, there is no need to introduce simulation and associated errors
% once we have fixed the shock process that generates randomness.
%
% The code here works when we are looking for the distribution of f(a,z),
% where a'(a,z), meaning that the a next period is determined by _a_ last
% period and some shock. Given this, the _a'_ is fixed for all _z'_. If
% however, the outcome of interest is such that: y'(y,z,z'), meaning that
% y' is different depending on realized z', the code below does not work,
% rather, this code
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_ds.html
% ff_iwkz_ds> should be used.
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
% * $p(Y=y, a) = \sum_{z} \left( 1\left\{Y(a,z)=y\right\} \cdot p(a,z) \right)$
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
%    [param_map, support_map] = ffs_az_set_default_param(it_param_set);
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
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_ds_post_stats.html ff_az_ds_post_stats>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html fft_disc_rand_var_stats>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html fft_disc_rand_var_mass2outcomes>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2covar.html fft_disc_rand_var_mass2covar>
%
% @seealso
%
% * derive distribution f(y'(y,z)) one asset *loop*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html ff_az_ds>
% * derive distribution f(y'({x,y},z)) two assets *loop*: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_ds.html ff_akz_ds>
% * derive distribution f(y'({x,y},z, *z'*)) two assets *loop*: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_ds.html ff_iwkz_ds>
% * derive distribution f(y'({y},z)) or f(y'({x,y},z)) *vectorized*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vec.html ff_az_ds_vec>
% * derive distribution f(y'({y},z, *z'*)) or f(y'({x,y},z, *z'*)) *vectorized*: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_ds_vec.html ff_iwkz_ds_vec>
% * derive distribution f(y'({y},z)) or f(y'({x,y},z)) *semi-analytical*: <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vecsv.html ff_az_ds_vecsv>
% * derive distribution f(y'({y},z, *z'*)) or f(y'({x,y},z, *z'*)) *semi-analytical*: <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_ds_vecsv.html ff_iwkz_ds_vecsv>
%

%% Default
% Program can be externally invoked with _az_, _abz_ or various other
% programs. By default, program invokes using _az_ model programs:
%
% # it_subset = 5 is basic invoke quick test
% # it_subset = 6 is invoke full test
% # it_subset = 7 is profiling invoke
% # it_subset = 8 is matlab publish
% # it_subset = 9 is invoke operational (only final stats) and coh graph
%

if (~isempty(varargin))
    
    % if invoked from outside override fully
    [param_map, support_map, armt_map, func_map, result_map] = varargin{:};

else
    % default invoke
    close all;
    
    it_param_set = 8;
    bl_input_override = true;

    % 1. Generate Parameters
    [param_map, support_map] = ffs_az_set_default_param(it_param_set);

    % Note: param_map and support_map can be adjusted here or outside to override defaults
    % param_map('it_a_n') = 750;
    % param_map('it_z_n') = 15;

    % 2. Generate function and grids
    [armt_map, func_map] = ffs_az_get_funcgrid(param_map, support_map, bl_input_override); % 1 for override
    
    % 3. Solve value and policy function using az_vf_vecsv, if want to solve
    % other models, solve outside then provide result_map as input
    [result_map] = ff_az_vf_vecsv(param_map, support_map, armt_map, func_map);    
    
end

%% Parse Parameters

% append function name
st_func_name = 'ff_az_ds';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

% result_map
% ar_st_pol_names is from section _Process Optimal Choices_ in the value
% function code.
params_group = values(result_map, {'cl_mt_pol_a'});
[cl_mt_pol_a] = params_group{:};
mt_pol_a = deal(cl_mt_pol_a{1});

% armt_map
params_group = values(armt_map, {'ar_a', 'mt_z_trans'});
[ar_a, mt_z_trans] = params_group{:};

% param_map
params_group = values(param_map, {'it_a_n', 'it_z_n'});
[it_a_n, it_z_n] = params_group{:};
params_group = values(param_map, {'it_maxiter_dist', 'fl_tol_dist'});
[it_maxiter_dist, fl_tol_dist] = params_group{:};

% support_map
params_group = values(support_map, {'bl_profile_dist', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time', 'bl_display_dist', 'it_display_every'});
[bl_profile_dist, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time, bl_display_dist, it_display_every] = params_group{:};

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

mt_dist_az_init = ones(length(ar_a),it_z_n)/length(ar_a)/it_z_n;
mt_dist_az_cur = mt_dist_az_init;
mt_dist_az_zeros = zeros(length(ar_a),it_z_n);

%% *f(a,z)*: Initialize Convergence Conditions

bl_histiter_continue = true;
it_iter = 0;
ar_dist_diff_norm = zeros([it_maxiter_dist, 1]);
mt_dist_perc_change = zeros([it_maxiter_dist, it_z_n]);

%% *f(a,z)*: Derive Stationary Distribution
% Iterate over the discrete joint random variable variables (a,z)
%
% We are looking for the distribution of: $p(a,z)$ where $a'(a,z)$, meaning
% that the a next period is determined by a last period and some shock.
% Given this, the $a'$ is fixed for all $z'$
%
% To make the code work for life-cycle model:
% # _mt_dist_az_init_: Initialize with potentially exogenous initial asset
% distribution
% # _mt_dist_az_: change mt_dist_az to tensor with a third dimension for
% age
% # at the beginning of the third loop over ar_z, get mass at current age,
% meaning: fl_cur_za_prob = ts_dist_az(it_a_prime_idx, it_zp_q, age)
% # at the end of the third loop over ar_z, add accumulated mass to next
% period, meaning: ts_dist_az(it_a_prime_idx, it_zp_q, age+1) =+ fl_zfromza
%

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
    for it_z_i = 1:it_z_n
        
        % loop 2: over endogenous states
        for it_a_j = 1:length(ar_a)
            
            % f(a'|a) = 1 for only one a'
            % in dynamic programming problem, had a loop over choices, now
            % already have optimal choices, do not need to loop
            fl_aprime = mt_pol_a(it_a_j, it_z_i);           
            it_a_prime_idx = find(ar_a == fl_aprime);
            
            % loop 3: loop over future shocks
            % E_{a,z}(f(a',z'|a,z)*f(a,z))
            for it_zp_q = 1:it_z_n
                
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
% function of the random variables (a,z), using f(a,z). We call function
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_ds_post_stats.html
% ff_az_ds_post_stats> which uses
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html
% fft_disc_rand_var_stats> and
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html
% fft_disc_rand_var_mass2outcomes> to compute various statistics of
% interest.

bl_input_override = true;
result_map = ff_az_ds_post_stats(support_map, result_map, mt_dist_az, bl_input_override);

end
