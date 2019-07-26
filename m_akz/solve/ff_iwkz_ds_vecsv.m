%% Derive Distributions for Risky + Safe Asets + Interpolated Distribution (Analytical)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [result_map] = ff_iwkz_ds_vecsv(varargin)
%% FF_IWKZ_DS finds the stationary asset distributions
% Building on the Two Assets Two-Step Interpolated Dynamic Programming
% Problem
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_iwkz_vf_vecsv.html
% ff_iwkz_vf_vecsv>, here we solve for the asset distribution. This version
% of the program is semi-analytical.
%
% This is the two-stage with interpolation version of
% <https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_ds_vecsv.html
% ff_akz_ds_vecsv>. See that file for additional descriptions and
% comparisons. These two functions are nearly identical
%
% The code here works when we are looking for the distribution of f(a,z),
% where a'(a,z,z'), meaning that the a next period is determined by a last
% period and some shock last period as well as shock this period. a here is
% cash-on-hand. This contrasts with
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vecsv.html
% ff_az_ds_vecsv>, which works for a'(a,z), a' can not be a function of z'.
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
%    ff_iwkz_ds_vecsv(param_map, support_map);
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

    it_param_set = 6;
    st_akz_or_iwkz = 'iwkz';

    % 1. Generate Parameters
    [param_map, support_map] = ffs_akz_set_default_param(it_param_set);

    % Note: param_map and support_map can be adjusted here or outside to override defaults
    % param_map('it_w_n') = 50;
    % param_map('it_z_n') = 15;
    
    param_map('fl_beta') = 0.90684;

    param_map('st_analytical_stationary_type') = 'eigenvector';

    % 2. Generate function and grids
    [armt_map, func_map] = ffs_akz_get_funcgrid(param_map, support_map); % 1 for override

    % 3. Solve value and policy function using ff_iwkz_vf_vecsv
    if (strcmp(st_akz_or_iwkz, 'iwkz'))
        [result_map] = ff_iwkz_vf_vecsv(param_map, support_map, armt_map, func_map);
    end
end

%% Parse Parameters

% append function name
st_func_name = 'ff_iwkz_ds_vecsv';
support_map('st_profile_name_main') = [st_func_name support_map('st_profile_name_main')];
support_map('st_mat_name_main') = [st_func_name support_map('st_mat_name_main')];
support_map('st_img_name_main') = [st_func_name support_map('st_img_name_main')];

% result_map
% ar_st_pol_names is from section _Process Optimal Choices_ in the value
% function code.
params_group = values(result_map, {'cl_mt_pol_a', 'cl_mt_pol_k'});
[cl_mt_pol_a, cl_mt_pol_k] = params_group{:};
[mt_pol_a, mt_pol_k] = deal(cl_mt_pol_a{1}, cl_mt_pol_k{1});

% Get Model Name
params_group = values(param_map, {'st_model'});
[st_model] = params_group{:};
% param_map
params_group = values(param_map, {'it_z_n'});
[it_z_n] = params_group{:};

% func_map
params_group = values(func_map, {'f_coh'});
[f_coh] = params_group{:};

% armt_map
params_group = values(armt_map, {'mt_z_trans', 'ar_interp_coh_grid'});
[mt_z_trans, ar_interp_coh_grid] = params_group{:};
if (ismember(st_model, ["ipwkbzr"]))
    params_group = values(armt_map, {'ar_z_r_borr_mesh_wage_w1r2', 'ar_z_wage_mesh_r_borr_w1r2'});
    [ar_z_r_borr_mesh_wage_w1r2, ar_z_wage_mesh_r_borr_w1r2] = params_group{:};
    params_group = values(param_map, {'it_z_wage_n', 'fl_z_r_borr_n'});
    [it_z_wage_n, fl_z_r_borr_n] = params_group{:};
elseif (ismember(st_model, ["ipwkbzr_fibs"]))
    params_group = values(armt_map, {'ar_z_r_infbr_mesh_wage_w1r2', 'ar_z_wage_mesh_r_infbr_w1r2'});
    [ar_z_r_borr_mesh_wage_w1r2, ar_z_wage_mesh_r_borr_w1r2] = params_group{:};
    params_group = values(param_map, {'it_z_wage_n', 'fl_z_r_infbr_n'});
    [it_z_wage_n, fl_z_r_borr_n] = params_group{:};
else
    params_group = values(armt_map, {'ar_z'});
    [ar_z] = params_group{:};
end

% param_map
params_group = values(param_map, {'st_analytical_stationary_type'});
[st_analytical_stationary_type] = params_group{:};

% support_map
params_group = values(support_map, {'bl_profile_dist', 'st_profile_path', ...
    'st_profile_prefix', 'st_profile_name_main', 'st_profile_suffix',...
    'bl_time'});
[bl_profile_dist, st_profile_path, ...
    st_profile_prefix, st_profile_name_main, st_profile_suffix, ...
    bl_time] = params_group{:};

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
it_exostates_n = it_z_n;

%% B. Solve for Index
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
if (ismember(st_model, ["ipwkbzr"]))
    mt_coh_prime = f_coh(ar_z_r_borr_mesh_wage_w1r2, ar_z_wage_mesh_r_borr_w1r2, ...
                        mt_pol_a(:), mt_pol_k(:));
elseif (ismember(st_model, ["ipwkbzr_fibs"]))
    % mt_pol_a includes interest rates
    mt_coh_prime = f_coh(ar_z_wage_mesh_r_borr_w1r2, mt_pol_a(:), mt_pol_k(:));
else
    mt_coh_prime = f_coh(ar_z, mt_pol_a(:), mt_pol_k(:));
end


% 2. *mt_coh_prime_on_grid_idx* is (coh_n x z_n) by (z_n):
% index for coh'(a,k,z')
[~, ar_coh_prime_on_grid_idx] = min(abs(mt_coh_prime(:)' - ar_interp_coh_grid'));
mt_coh_prime_on_grid_idx = reshape(ar_coh_prime_on_grid_idx, size(mt_coh_prime));

%% C. Expand Index so Matches Full States Index Dimension
% The index above matches the index in the cash-on-hand grid, but now, the
% state space is cash-on-hand jointly with shocks, that is the full states
% markov's states. So if there are two shocks and two cash-on-hand grid
% points, the cash-on-hand grid points would have been [1,2] and [1,2], but
% depending on which z' they match up to, they would now be [1,2] if
% matching to the first z', and [3,4] if matching to the second z'.

% mt_pol_idx_mesh_max is (NxM) by M, mt_pol_idx is N by M
mt_pol_idx_mesh_max = mt_coh_prime_on_grid_idx + (0:1:(it_exostates_n-1))*it_endostates_n;

%% D. Transition Probabilities from (M by M) to (NxM) by M
% Probability comes from the shock transition matrix, which is now
% duplicated for all cash-on-hand grid elements

mt_trans_prob = reshape(repmat(mt_z_trans(:)', ...
    [it_endostates_n, 1]), [it_endostates_n*it_exostates_n, it_exostates_n]);

%% E. Fill mt_pol_idx_mesh_idx to mt_full_trans_mat SPARSE
% Try to always use sparse matrix, unless grid sizes very small, keeping
% non-sparse code here for comparison. Sparse matrix is important for
% allowing the code to be fast and memory efficient. Otherwise this method
% is much slower than iterative method.

i = mt_pol_idx_mesh_max(:);
j = repmat((1:1:it_endostates_n*it_exostates_n),[1,it_exostates_n])';
v = mt_trans_prob(:);
m = it_endostates_n*it_exostates_n;
n = it_endostates_n*it_exostates_n;
mt_full_trans_mat = sparse(i, j, v, m, n);

%% F. Stationary Distribution *Method A*, Eigenvector Approach
% Given that markov chain we have constructured for all state-space
% elements, we can now find the stationary distribution using standard
% <https://en.wikipedia.org/wiki/Markov_chain#Stationary_distribution_relation_to_eigenvectors_and_simplices
% eigenvector> approach. See
% <https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_ds_vecsv.html
% ff_az_ds_vecsv> for additional methods using the full states markov
% structure.

if (strcmp(st_analytical_stationary_type, 'eigenvector'))
    [V, ~] = eigs(mt_full_trans_mat,1,1);
    ar_stationary = V/sum(V);
end

%% G. Stationary Vector to Stationary Matrix in Original Dimensions

mt_dist_akz = reshape(ar_stationary, size(mt_pol_a));


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
