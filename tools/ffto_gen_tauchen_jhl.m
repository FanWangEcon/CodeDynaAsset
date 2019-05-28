%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is modified based on
% https://sites.google.com/site/janhanneslang/programs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [s, Pi, stationary, z_vector] = fo_st_tauchen_jhl(mu ,rho, sig, N)

m       = 3;
s       = zeros(N,1);
Pi      = zeros(N,N);
s(1)    = mu/(1-rho) - m*sqrt(sig^2/(1-rho^2));
s(N)    = mu/(1-rho) + m*sqrt(sig^2/(1-rho^2));
step    = (s(N)-s(1))/(N-1);

for i=2:(N-1)
   s(i) = s(i-1) + step;
end

for j = 1:N
    for k = 1:N
        if k == 1
            Pi(j,k) = cdf_normal((s(1) - mu - rho*s(j) + step/2) / sig);
        elseif k == N
            Pi(j,k) = 1 - cdf_normal((s(N) - mu - rho*s(j) - step/2) / sig);
        else
            Pi(j,k) = cdf_normal((s(k) - mu - rho*s(j) + step/2) / sig) - ...
                      cdf_normal((s(k) - mu - rho*s(j) - step/2) / sig);
        end
    end
end

stationary = Pi^1000;
stationary = stationary(1,:);
Lagg = stationary*exp(s);
z_vector = exp(s')/Lagg;

end

function c = cdf_normal(x)
    c = 0.5 * erfc(-x/sqrt(2));
end
