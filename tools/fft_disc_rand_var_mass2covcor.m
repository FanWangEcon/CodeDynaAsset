%% Compute Covariance, Correlation for cov(x,y) given X(a,z), Y(a,z) and f(a,z)
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [fl_cov_xy, fl_cor_xy] = fft_disc_rand_var_mass2covcor(varargin)
%% FFT_DISC_RAND_VAR_MASS2COVCOR find cov(x,y) given X(a,z), Y(a,z) and f(a,z)
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
%
% @param st_var_name string name of the variable (choice/outcome) been analyzed
%
% @param mt_choice_bystates matrix N by M of choices along two dimensions,
% N could be endogenous states, M could be exogenous shocks, or vice-versa
%
% @param mt_dist_bystates matrix N by M of probability mass on states, N
% could be endogenous states, M could be exogenous shocks, or vice versa
%
% @return tb_choice_drv_cur_byY table table containing two columns, unique
% outcomes/choices y from y(a,z) and probability mass associated with each
% y f(y)
%
% @return ar_choice_prob_byY table array probability mass associated with each
% y f(y), second column from tb_choice_drv_cur_byY, dimension unknown,
% determined by y(a,z) function
%
% @return ar_choice_unique_sorted_byY table array unique Ys, dimension
% unknown, determined by y(a,z) function
%
% @return mt_choice_prob_byYZ matrix f(y,z), meaning for y outcomes along
% the column dimension.
%
% @return mt_choice_prob_byYA matrix f(y,a), meaning for y outcomes along
% the row dimension.
%
% @seealso
%
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_stats.html fft_disc_rand_var_stats>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2outcomes.html fft_disc_rand_var_mass2outcomes>
% * <https://fanwangecon.github.io/CodeDynaAsset/tools/html/fft_disc_rand_var_mass2covcor.html fft_disc_rand_var_mass2covcor>
%

%% Default
% use binomial as test case, z maps to binomial win prob, remember binom
% approximates normal.

if (~isempty(varargin))
    
    % if invoked from outside overrid fully
    [covvar_input_map] = varargin{:};
    bl_display_drvm2covcor = false;
    
else
    
    clear all;
    close all;
    
    it_states = 6;
    it_shocks = 5;
    fl_binom_n = it_states-1;
    ar_binom_p = (1:(it_shocks))./(it_shocks+2);
    ar_binom_x = (0:1:(it_states-1)) - 3;
    
    % f(z)
    ar_binom_p_prob = binopdf(0:(it_shocks-1), it_shocks-1, 0.5);
    % f(a,z), mass for a, z
    mt_dist_bystates = zeros([it_states, it_shocks]);
    for it_z=1:it_shocks
        % f(a|z)
        f_a_condi_z = binopdf(ar_binom_x - min(ar_binom_x), fl_binom_n, ar_binom_p(it_z));
        % f(z)
        f_z = ar_binom_p_prob(it_z);
        % f(a,z)=f(a|z)*f(z)
        mt_dist_bystates(:, it_z) = f_a_condi_z*f_z;
    end
    
    % x(a,z), some non-smooth structure
    rng(123);
    mt_choice_x_bystates = ar_binom_x' - 0.01*ar_binom_x'.^2  + ar_binom_p - 0.5*ar_binom_p.^2 + rand([it_states, it_shocks]);
    mt_choice_x_bystates = round(mt_choice_x_bystates*3);

    % y(a,z), some non-smooth structure
    rng(456);
    mt_choice_y_bystates = 10 -(mt_choice_x_bystates) + 15*(rand([it_states, it_shocks])-0.5);
        
    % Obtain mean and sd
    st_cur_output_key = 'x_outcome';
    [ar_choice_prob_byX, ar_choice_unique_sorted_byX, ~, ~] = ...
        fft_disc_rand_var_mass2outcomes(st_cur_output_key, mt_choice_x_bystates, mt_dist_bystates);
    [ds_stats_x_map] = fft_disc_rand_var_stats(st_cur_output_key, ar_choice_unique_sorted_byX', ar_choice_prob_byX');
    fl_choice_x_mean = ds_stats_x_map('fl_choice_mean');
    fl_choice_x_sd = ds_stats_x_map('fl_choice_sd');
    
    st_cur_output_key = 'y_outcome';
    [ar_choice_prob_byY, ar_choice_unique_sorted_byY, ~, ~] = ...
        fft_disc_rand_var_mass2outcomes(st_cur_output_key, mt_choice_y_bystates, mt_dist_bystates);
    [ds_stats_y_map] = fft_disc_rand_var_stats(st_cur_output_key, ar_choice_unique_sorted_byY', ar_choice_prob_byY');
    fl_choice_y_mean = ds_stats_y_map('fl_choice_mean');
    fl_choice_y_sd = ds_stats_y_map('fl_choice_sd');
        
    % display
    bl_display_drvm2covcor = true;
    
    % Collect
    covvar_input_map = containers.Map('KeyType','char', 'ValueType','any');
    covvar_input_map('mt_choice_x_bystates') = mt_choice_x_bystates;
    covvar_input_map('mt_choice_y_bystates') = mt_choice_y_bystates;
    covvar_input_map('mt_dist_bystates') = mt_dist_bystates;
    covvar_input_map('fl_choice_x_mean') = fl_choice_x_mean;
    covvar_input_map('fl_choice_x_sd') = fl_choice_x_sd;
    covvar_input_map('fl_choice_y_mean') = fl_choice_y_mean;
    covvar_input_map('fl_choice_y_sd') = fl_choice_y_sd;
end

%% Parse Parameters

% probability over a and z
params_group = values(covvar_input_map, {'mt_dist_bystates'});
[mt_dist_bystates] = params_group{:};

% x and y outcomes
params_group = values(covvar_input_map, {'mt_choice_x_bystates', 'mt_choice_y_bystates'});
[mt_choice_x_bystates, mt_choice_y_bystates] = params_group{:};

% x and y stats
params_group = values(covvar_input_map, {'fl_choice_x_mean', 'fl_choice_x_sd', ...
                                         'fl_choice_y_mean', 'fl_choice_y_sd'});
[fl_choice_x_mean, fl_choice_x_sd, fl_choice_y_mean, fl_choice_y_sd] = params_group{:};

%% 1. Compute Covariance

mt_x_devi_from_mean = (mt_choice_x_bystates - fl_choice_x_mean);
mt_y_devi_from_mean = (mt_choice_y_bystates - fl_choice_y_mean);
mt_x_y_multiply = (mt_x_devi_from_mean).*(mt_y_devi_from_mean);
mt_cov_component_weighted = mt_dist_bystates.*(mt_x_y_multiply);
fl_cov_xy = sum(mt_cov_component_weighted, 'all');

%% 2. Compute Correlation

fl_cor_xy = fl_cov_xy/(fl_choice_x_sd*fl_choice_y_sd);

%% Display
if (bl_display_drvm2covcor)
        
    fft_container_map_display(covvar_input_map, 25, 15);
        
    covvar_output_map = containers.Map('KeyType','char', 'ValueType','any');
    covvar_output_map('mt_x_devi_from_mean') = mt_x_devi_from_mean;
    covvar_output_map('mt_y_devi_from_mean') = mt_y_devi_from_mean;
    covvar_output_map('mt_x_y_multiply') = mt_x_y_multiply;
    covvar_output_map('mt_cov_component_weighted') = mt_cov_component_weighted;
    
    fft_container_map_display(covvar_output_map, 25, 15);    
    
    disp('fl_cov');
    disp(fl_cov_xy);
    
    disp('fl_cor');
    disp(fl_cor_xy);
    
end
end