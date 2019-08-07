%% Generate Discrete Random Variable
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [ar_drv_ele, ar_drv_prb] = fft_gen_discrete_var(varargin)
%% FFT_GEN_DISCRETE_VAR generate discrete random variables
% could be interest rates. Given parameters, generate discrete random
% variable, translating to other level/units elsewhere potentially.
%
% @st_drv_ele_type string scalar values associated with element of the
% sample space of the discrete random variable, possibilities are:
% 
% # unif: uniform between _fl_max_ and _fl_min_, equi-distance. 
% # seg3: between _fl_max_ and _fl_min_, three gaps, increasing for each
% segment, might not reach _fl_max_
% # logspace: logspace function between _fl_max_ and _fl_min_
%
% @st_drv_prb_type string distributional functions, possibilities are:
%
% # poiss: top truncated poisson distribution (non-symmetric, approximates exponential)
% # binom: binomial (can be symmetric, approximates normal) 
% # unif: uniform distribution
%
% @params_map container an container of parameters if the parameter is not covered by the
% explicitly named parameters
%
% @return ar_drv_ele array discrete random variable outcomes
%
% @return ar_drv_prb array discrete random variable probabilities
%

%% Parameters Defaults

bl_input_override = 0;
if (length(varargin) == 2)
    bl_input_override = varargin{2};
end

if (bl_input_override)
    
    [param_dsv_map, ~] = varargin{:};
    
    bl_display_dvg = false;
    
else
    
    param_dsv_map = containers.Map('KeyType','char', 'ValueType','any');
    
    param_dsv_map('fl_binom_p') = 0.25;
    param_dsv_map('fl_poiss_mean') = 5;
    param_dsv_map('fl_logspace_adj') = 0.075;
    
    param_dsv_map('st_drv_ele_type') = 'unif';
    param_dsv_map('st_drv_prb_type') = 'poiss';
    param_dsv_map('fl_max') = 0.095;
    param_dsv_map('fl_min') = 0.025;
    param_dsv_map('fl_n') = 3;
    
    bl_display_dvg = true;
end

%% Retrieve Parameters from Map

params_group = values(param_dsv_map, {'st_drv_ele_type', 'st_drv_prb_type', 'fl_max', 'fl_min', 'fl_n'});
[st_drv_ele_type, st_drv_prb_type, fl_max, fl_min, fl_n] = params_group{:};

%% Gereate Interest Rates Shocks
if (fl_max == fl_min)
%% Single Point Shock

    ar_drv_ele = linspace(fl_min, fl_max, 1);
    ar_drv_prb = zeros(size(ar_drv_ele)) + 1;

else
    
%% Multiple Points Shock    
    %% Discrete Outcomes: Uniform
    fl_gap = (fl_max - fl_min)/(fl_n-1);
    
    if (strcmp(st_drv_ele_type, 'unif'))
        ar_drv_ele = fl_min:fl_gap:fl_max;
    end
    
    %% Discrete Outcomes: three segments of increasing gaps
    
    if (strcmp(st_drv_ele_type, 'seg3'))
        
        fl_most_least_gap = (fl_max - fl_min);
        fl_most_least_seg3_interval = fl_most_least_gap/3;
        
        ar_seg1 = fl_min:fl_gap:fl_most_least_seg3_interval;
        ar_seg2 = max(ar_seg1):(fl_gap*2):fl_most_least_seg3_interval*2;
        ar_seg3 = max(ar_seg2):(fl_gap*3):fl_max;
        
        ar_drv_ele =[ar_seg1, ar_seg2, ar_seg3];
    end
    
    %% Discrete Outcomes: log space
    
    if (strcmp(st_drv_ele_type, 'logspace'))
        fl_logspace_adj = param_dsv_map('fl_logspace_adj');
        ar_drv_ele = logspace(log10(fl_min + fl_logspace_adj), ...
            log10(fl_max + fl_logspace_adj), fl_n) - fl_logspace_adj;
    end
    
    %% Probability: Uniform
    
    if (strcmp(st_drv_prb_type, 'unif'))
        ar_drv_prb = zeros(size(ar_drv_ele)) + 1/length(ar_drv_ele);
    end
    
    %% Probability: binom
    
    if (strcmp(st_drv_prb_type, 'binom'))
        fl_binom_p = param_dsv_map('fl_binom_p') ;
        it_poiss_n = length(ar_drv_ele)-1;
        ar_poiss_x = 0:1:it_poiss_n;
        ar_drv_prb = binopdf(ar_poiss_x, it_poiss_n, fl_binom_p);
    end
    
    %% Probability: poisson
    % divide by sum because poisson max n is infinity
    
    if (strcmp(st_drv_prb_type, 'poiss'))
        fl_poiss_mean = param_dsv_map('fl_poiss_mean');
        it_poiss_n = length(ar_drv_ele)-1;
        ar_poiss_x = 0:1:it_poiss_n;
        ar_drv_prb = poisspdf(ar_poiss_x, fl_poiss_mean);
        ar_drv_prb = ar_drv_prb./sum(ar_drv_prb);
    end
    
end

%% Display

if (bl_display_dvg)
    disp('ar_drv_ele');
    disp(ar_drv_ele);
    disp('ar_drv_prb');
    disp(ar_drv_prb);
end

end
