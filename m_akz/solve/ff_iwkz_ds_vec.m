%% Derive Distributions for Risky + Safe Asets + Interpolated Distribution (Vectorized)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_iwkz_ds_vec(varargin)
%% FF_IWKZ_DS_VEC finds the stationary asset distributions
% Building on the Two Assets Two-Step Interpolated Dynamic Programming
% Problem
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_vf_vecsv.html
% ff_iwkz_vf_vecsv>, here we solve for the asset distribution. This version
% of the program is vectorized
%
% This is the two-stage with interpolation version of
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_ds_vec.html
% ff_akz_ds_vec>. See that file for additional descriptions and
% comparisons. These two functions are nearly identical
%
% The code here works when we are looking for the distribution of f(a,z),
% where a'(a,z,z'), meaning that the a next period is determined by a last
% period and some shock last period as well as shock this period. a here is
% cash-on-hand. This contrasts with
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html
% ff_az_ds>, which works for a'(a,z), a' can not be a function of z'.
%
% @example
%
%    % Get Default Parameters
%    it_param_set = 6;
%    [param_map, support_map] = ffs_az_set_default_param(it_param_set);
%    % Change Keys in param_map
%    param_map('it_w_n') = 750;
%    param_map('it_ak_n') = param_map('it_w_n');
%    param_map('it_z_n') = 11;
%    param_map('fl_a_max') = 100;
%    param_map('fl_w') = 1.3;
%    % Change Keys support_map
%    support_map('bl_display') = false;
%    support_map('bl_post') = true;
%    support_map('bl_display_final') = false;
%    % Call Program with external parameters that override defaults
%    ff_iwkz_ds_vec(param_map, support_map);
%
% @include
%
% * <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_vf_vecsv.html ff_wkz_vf_vecsv>
% * <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_ds_post_stats.html ff_az_ds_post_stats>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html fft_disc_rand_var_stats>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html fft_disc_rand_var_mass2outcomes>
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
    st_akz_or_iwkz = 'iwkz';

    % 1. Generate Parameters
    [param_map, support_map] = ffs_akz_set_default_param(it_param_set);

    % Note: param_map and support_map can be adjusted here or outside to override defaults
    % param_map('it_w_n') = 50;
    % param_map('it_z_n') = 15;

    % 2. Generate function and grids
    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map); % 1 for override

    % 3. Solve value and policy function using ff_iwkz_vf_vecsv
    if (strcmp(st_akz_or_iwkz, 'iwkz'))
        [result_map] = ff_iwkz_vf_vecsv(param_map, support_map, armt_map, func_map);
    end
end

%% Parse Parameters

% append function name
st_func_name = 'ff_iwkz_ds_vec';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

% result_map
% ar_st_pol_names is from section _Process Optimal Choices_ in the value
% function code.
params_group = values(result_map, {'cl_mt_pol_a', 'cl_mt_pol_k'});
[cl_mt_pol_a, cl_mt_pol_k] = params_group{:};
[mt_pol_a, mt_pol_k] = deal(cl_mt_pol_a{1}, cl_mt_pol_k{1});

% func_map
params_group = values(func_map, {'f_coh'});
[f_coh] = params_group{:};

% armt_map
params_group = values(armt_map, {'mt_z_trans', 'ar_z'});
[mt_z_trans, ar_z] = params_group{:};
params_group = values(armt_map, {'ar_interp_coh_grid'});
[ar_interp_coh_grid] = params_group{:};

% param_map
params_group = values(param_map, {'it_z_n', 'it_maxiter_dist', 'fl_tol_dist'});
[it_z_n, it_maxiter_dist, fl_tol_dist] = params_group{:};

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

%% A. Get Size of Endogenous and Exogenous State

it_endostates_n = length(ar_interp_coh_grid);
it_exostates_n = length(ar_z);

%% B. Initialize Output Matrixes
% Initialize the distribution to be uniform

mt_dist_akz_init = ones(it_endostates_n,it_exostates_n)/it_endostates_n/it_exostates_n;
mt_dist_akz_cur = mt_dist_akz_init;
mt_dist_akz_zeros = zeros(it_endostates_n,it_exostates_n);

%% C. Initialize Convergence Conditions

bl_histiter_continue = true;
it_iter = 0;
ar_dist_diff_norm = zeros([it_maxiter_dist, 1]);
mt_dist_perc_change = zeros([it_maxiter_dist, it_z_n]);

%% D. Solve for Index
% The model is solved by interpolating over cash-on-hand. The optimal
% choices do not map to specific points on the cash-on-hand grid. Find the
% index of the cash-on-hand vector that is the closest to the
% coh'(a'(coh,z),k'(coh,z),z'). 
%
% Since we have *z_n* elements of shocks, and *coh_n* elements of the
% cash-on-hand grid, there are (coh_n x z_n) possible combinations of
% states at period t. In period t+1, there are (coh_n x z_n) by (z_n)
% possible/reachable cash-on-hand points. We find the index of all these
% reachable coh' points on the interpolation cash-on-hand grid.
%

% 1. *mt_coh_prime* is (coh_n x z_n) by (z_n)
% coh'(z', a'(coh,z), k'(coh,z))    
mt_coh_prime = f_coh(ar_z, mt_pol_a(:), mt_pol_k(:));

% 2. *mt_coh_prime_on_grid_idx* is (coh_n x z_n) by (z_n):
% index for coh'(a,k,z')
mt_coh_prime_on_grid_idx = zeros(size(mt_coh_prime));
for it_zprime_ctr=1:size(mt_coh_prime, 2)
    ar_coh_prime = mt_coh_prime(:,it_zprime_ctr);
    [~, ar_coh_prime_on_grid_idx] = min(abs(ar_coh_prime(:)' - ar_interp_coh_grid'));
    mt_coh_prime_on_grid_idx(:,it_zprime_ctr) = ar_coh_prime_on_grid_idx;
end

%% E. Solve for Unique Index
% For each z', there are (coh_n x z_n) possible coh'(z') reachable points,
% which have been converted to index: *mt_coh_prime_on_grid_idx* along the
% cash-on-hand grid in the previous code segment. Now, we find the number
% of unique coh' grid points among the (coh_n x z_n) grid indexes:
% *ar_idx_of_unique*. We also find the positions of these unique indexes in
% the full (coh_n x z_n) grid: *ar_idx_full*
% 

cl_ar_idx_full = cell([it_exostates_n, 1]);
cl_ar_idx_of_unique = cell([it_exostates_n, 1]);

for it_z_i = 1:it_exostates_n
    
    % 5. Cumulative probability received at state from zi
    [ar_idx_full, ~, ar_idx_of_unique] = unique(mt_coh_prime_on_grid_idx(:, it_z_i));
    
    cl_ar_idx_full{it_z_i} = ar_idx_full;
    cl_ar_idx_of_unique{it_z_i} = ar_idx_of_unique;
    
end

%% F. Derive Stationary Distribution
% Iterate until convergence

while (bl_histiter_continue)

    it_iter = it_iter + 1;

    %% F1. Iterate over z' Shocks
    % The code below loops over future states, note that the structure here is
    % significant different from
    % <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_ds.html
    % ff_akz_ds>, where we looped over current shocks. Here we loop over future
    % shocks because coh' is a function of z' as well as choices last period. 
    %
    
    % 1. initialize empty
    mt_dist_akz = mt_dist_akz_zeros;
            
    % 2. loop over next period exo shocks
    for it_z_i = 1:it_exostates_n

        % 3. *ar_zi_prob* is (coh_n x z_n) by (1):
        % mt_z_trans(:, it_z_i)': for all z today, prob(z'|z) fixing z' 
        % overall: f(coh,z)*f(z'|z) fixing z' for all z
        mt_zi_prob = mt_dist_akz_cur .* mt_z_trans(:,it_z_i)';
        ar_zi_prob = mt_zi_prob(:);
        
        % 4. Cumulative probability received at state from zi
        mt_zi_cumu_prob = accumarray(cl_ar_idx_of_unique{it_z_i}, ar_zi_prob);

        % 5. Adding up
        mt_dist_akz(cl_ar_idx_full{it_z_i}, it_z_i) = mt_zi_cumu_prob +  mt_dist_akz(cl_ar_idx_full{it_z_i}, it_z_i);
    end

    %% F2. Check Tolerance and Continuation

    % Difference across iterations
    ar_dist_diff_norm(it_iter) = norm(mt_dist_akz - mt_dist_akz_cur);
    mt_dist_perc_change(it_iter, :) = sum((mt_dist_akz ~= mt_dist_akz))/it_endostates_n;

    % Update
    mt_dist_akz_cur = mt_dist_akz;

    % Print Iteration Results
    if (bl_display_dist && (rem(it_iter, it_display_every)==0))
        fprintf('Dist it_iter:%d, fl_dist_diff:%d\n', it_iter, ar_dist_diff_norm(it_iter));
        tb_hist_iter = array2table([sum(mt_dist_akz_cur,1); std(mt_dist_akz_cur,1); ...
                                    mt_dist_akz_cur(1,:); mt_dist_akz_cur(it_endostates_n,:)]);
        tb_hist_iter.Properties.VariableNames = strcat('z', string((1:size(mt_dist_akz,2))));
        tb_hist_iter.Properties.RowNames = {'mdist','sddist', 'Ldist', 'Hdist'};
        disp('mdist = sum(mt_dist_akz_cur,1) = sum_{a,k}(p({a,k})|z)')
        disp('sddist = std(mt_pol_ak_cur,1) = std_{a,k}(p({a,k})|z)')
        disp('Ldist = mt_dist_akz_cur(1,:) = p(min({a,k})|z)')
        disp('Hdist = mt_dist_akz_cur(it_a_n,:) = p(max({a,k})|z)')
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


%% *f(y), f(c), f(a), f(k)*: Generate Key Distributional Statistics for Each outcome
% Having derived f({a,k},z) the probability mass function of the joint discrete
% random variables, we now obtain distributional statistics. Note that we
% know f({a,k},z), and we also know relevant policy functions a'(a,k,z), k'(a,k,z),
% or other policy functions. We can simulate any choices that are a
% function of the random variables (coh(a,k),z), using f(coh(a,k),z). We call function
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_ds_post_stats.html
% ff_az_ds_post_stats> which uses
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html
% fft_disc_rand_var_stats> and
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html
% fft_disc_rand_var_mass2outcomes> to compute various statistics of
% interest.
result_map('mt_dist') = mt_dist_akz;
result_map = ff_az_ds_post_stats(support_map, result_map, mt_dist_akz);

end
