%% Add Zero to Array
% *back to <https://fanwangecon.github.io Fan>'s
% <https://fanwangecon.github.io/CodeDynaAsset/ Dynamic Assets Repository>
% Table of Content.*

%%
function [ar_a_wth_zero] = fft_array_add_zero(varargin)
%% FFT_ARRAY_ADD_ZERO Add zero to array without changing size
%

%% Default
% use binomial as test case, z maps to binomial win prob, remember binom
% approximates normal.

params_len = length(varargin);
bl_input_override = 0;
if (params_len == 2)
    bl_input_override = varargin{2};
end

if (bl_input_override)

    % if invoked from outside overrid fully
    [ar_a, ~] = varargin{:};
    [fl_a_min, fl_a_max, it_a_n] = deal(min(ar_a), max(ar_a), length(ar_a));
    bl_display_addzero = false;

else

    clear all;
    close all;
    
    [fl_a_min, fl_a_max, it_a_n] = deal(-20, 50, 11);
    ar_a = linspace(fl_a_min, fl_a_max, it_a_n);
    bl_display_addzero = true;

end

%% Add in Zero
if (ismember(0, ar_a))
    ar_a_wth_zero = sort(unique(ar_a));
else
    ar_a = linspace(fl_a_min, fl_a_max, it_a_n-1);
    
    if (ismember(0, ar_a))
        % add a mid point between 0 and lowest savings point above zero.
        ar_a_wth_zero = [(ar_a(2)-ar_a(1))/2, ar_a];
        ar_a_wth_zero = sort(unique(ar_a_wth_zero));
    else
        ar_a_wth_zero = [0 ar_a];
        ar_a_wth_zero = sort(unique(ar_a_wth_zero));        
    end
end

%% Display
if (bl_display_addzero)

    disp('ar_a');
    disp(ar_a);
    
    disp('ar_a_wth_zero');
    disp(ar_a_wth_zero);

end
end
