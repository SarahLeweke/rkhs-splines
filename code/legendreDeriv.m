function [legendreD, legendreP] = legendreDeriv(N, t)
% Leweke 2021

% The algorithm computes the frist derivative P_n'(t) of the Legendre Polynomials
% $P_n^\prime(t) for n=0,..,N

% returns legendreD (derivative) of size (N+1) x length(t)
% returns legendreP of size (N+1) x length(t)
    arguments
        N  (1,1) double {mustBeNumeric,mustBeReal,mustBeInteger} 
        t (:,1) double {mustBeNumeric,mustBeReal} 
    end

    % Linearise t (row vector)
    t = t';

    % Initialize Matrix
    legendreD = zeros(N+1, length(t)); 
    
    % Calculate Legendre Polynomials
    legendreP = legendrePoly(N,t);

    
    % algorithm can be found in for instance [Fengler Diss 2005 U Kaiserslautern Algo 15.3.2]
    % Test if norm(t) ~= 1
    idx = find(abs(abs(t)-1) > 1e-5);  
    if ~isempty(idx)
        for n = 1:N
            legendreD(n+1,idx)  = (n*t(idx).*legendreP(n+1,idx) - n*legendreP(n, idx))./(t(idx).^2 -1);
        end
    end

    % Test if t == 1
    idx = find(abs(t-ones(1,length(t))) <= 1e-5); 
    legendreD(1,idx) = zeros(1,length(t(idx)));
    legendreD(2,idx) = ones(1, length(t(idx)));
    if N > 2
        legendreD(3,idx) = 3*ones(1, length(t(idx)));
        for n = 3:N
           legendreD(n+1,idx) = 0.5*n*(n+1)*ones(1,length(t(idx)));
        end
    end
    
    % Test if t == -1
    idx = find(abs(t+ones(1,length(t))) <= 1e-5);   
    legendreD(1,idx) = zeros(1,length(t(idx)));
    legendreD(2,idx) = ones(1, length(t(idx)));
    if ~isempty(idx) && N > 2
        legendreD(3,idx) = -3*ones(1, length(t(idx)));
        for n = 3:N
           legendreD(n+1,idx) = (-1)^(n+1)*0.5*n*(n+1)*ones(1,length(t(idx)));
        end
    end
end