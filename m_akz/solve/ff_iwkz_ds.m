%% Derive Distributions for Risky + Safe Asets + Interpolated Distribution (Loop)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_iwkz_ds(varargin)
%% FF_IWKZ_DS finds the stationary asset distributions
% Building on the Two Assets Two-Step Interpolated Dynamic Programming
% Problem
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_vf_vecsv.html
% ff_iwkz_vf_vecsv>, here we solve for the asset distribution. This version
% of the program uses loops.
%
% This is the two-stage with interpolation version of
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_ds.html
% ff_akz_ds>. See that file, as well as
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html
% ff_az_ds> for additional descriptions and comparisons. These two
% functions are different. Specifically, the code here works when we are
% looking for the distribution of f(a,z), where a'(a,z,z'), meaning that
% the _a_ next period is determined by choices last period and some _z_
% last period, but also _z'_ this period. _a_ here specifically
% cash-on-hand. This contrasts with
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds.html
% ff_az_ds>, which works for a'(a,z), _a'_ can not be a function of _z'_.
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
%    ff_iwkz_ds(param_map, support_map);
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
    
    % override when called from outside
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
st_func_name = 'ff_iwkz_ds';
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

%% B. *f({a,k},z)*: Initialize Output Matrixes
% Initialize the distribution to be uniform

mt_dist_akz_init = ones(it_endostates_n,it_exostates_n)/it_endostates_n/it_exostates_n;
mt_dist_akz_cur = mt_dist_akz_init;
mt_dist_akz_zeros = zeros(it_endostates_n,it_exostates_n);

%% C. *f({a,k},z)*: Initialize Convergence Conditions

bl_histiter_continue = true;
it_iter = 0;
ar_dist_diff_norm = zeros([it_maxiter_dist, 1]);
mt_dist_perc_change = zeros([it_maxiter_dist, it_z_n]);

%% D. *f({a,k},z)*: Derive Stationary Distribution
% Iterate over the discrete joint random variable variables ({a,k},z)
%
% We are looking for the distribution of: $p({a,k},z)$ where $a'({a,k},z)$
% and $k'({a,k},z)$, meaning that $a'$ and $k'$ are determined by $a$ and
% $k$ and $z$, but not $z'$.
%

while (bl_histiter_continue)

    it_iter = it_iter + 1;

    %% *f({a,k},z)*: Iterate over Probability mass for Discrete Random Variable
    % compared to
    % <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_akz_vf.html
    % ff_akz_vf>, we basically have the same set of loops. There, there were
    % four loops, here there are three loops. We eliminated the loop over
    % next period choices, because here we already know optimal choices
    
    % initialize empty
    mt_dist_akz = mt_dist_akz_zeros;

    % loop 1: over exogenous states
    for it_z_i = 1:it_exostates_n

        % loop 2: over endogenous states
        for it_coh_grid_j = 1:it_endostates_n

            % f(a'|a) = 1 for only one a'
            % in dynamic programming problem, had a loop over choices, now
            % already have optimal choices, do not need to loop
            fl_aprime = mt_pol_a(it_coh_grid_j, it_z_i);
            fl_kprime = mt_pol_k(it_coh_grid_j, it_z_i);
            
            % loop 3: loop over future shocks
            % E_{coh,z}(f(coh'(a',k',z'),z'|coh,z)*f(coh,z))
            for it_zp_q = 1:it_exostates_n
                
                % A. Get the index that the index for coh' based on coh,z,z'
                % current shock
                fl_zprime = ar_z(it_zp_q);                
                % cash-on-hand next period which is a function also of z'
                fl_coh_prime = f_coh(fl_zprime, fl_aprime, fl_kprime);
                % next period index
                [~, it_coh_prime_on_grid_idx] = min(abs(fl_coh_prime - ar_interp_coh_grid));
                
                % B. prob of going to coh'(coh,z,z')
                % current probablity at (a,z)
                fl_cur_zak_prob = mt_dist_akz_cur(it_coh_grid_j, it_z_i);
                % f(z'|z) transition
                fl_ztoz_trans =  mt_z_trans(it_z_i, it_zp_q);
                % f(a',z'|a,z)*f({a,k},z)
                fl_zfromzak = fl_cur_zak_prob*fl_ztoz_trans;

                % cumulating
                mt_dist_akz(it_coh_prime_on_grid_idx, it_zp_q) = ...
                    mt_dist_akz(it_coh_prime_on_grid_idx, it_zp_q) + fl_zfromzak;
            end

        end

    end

    %% *f({a,k},z)*: Check Tolerance and Continuation

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
% know f({a,k},z), and we also know relevant policy functions a'(a,z), c(a,z),
% or other policy functions. We can simulate any choices that are a
% function of the random variables (a,z), using f({a,k},z). We call function
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_ds_post_stats.html
% ff_az_ds_post_stats> which uses
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html
% fft_disc_rand_var_stats> and
% <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html
% fft_disc_rand_var_mass2outcomes> to compute various statistics of
% interest.

result_map = ff_az_ds_post_stats(support_map, result_map, mt_dist_akz);

end
