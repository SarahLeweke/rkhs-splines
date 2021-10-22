function legendrePol = legendrePoly(N,t)
% Leweke 2021
% The algorithm computes the Legendre Polynomials
% P_n(t) for n=0,..,N

% returns legendrePol of size (N+1) x length(t)
    arguments
        N  (1,1) double {mustBeNumeric,mustBeReal,mustBeInteger} 
        t (:,1) double {mustBeNumeric,mustBeReal} 
    end

    % Linearize t (row vector)
    t = t';

    % Initialize Matrix
    legendrePol = zeros(N+1, length(t)); 

    % Initial values P_0(t) = 1, P_1(t) = t
    legendrePol(1,:) = ones(size(legendrePol,2),1);
    legendrePol(2,:) = t;

    
    % algorithm can be found in for instance [Fengler Diss 2005 U Kaiserslautern Algo 15.3.1]
    for n = 2:N
       legendrePol(n+1,:)  = ((2*n-1).*t.*legendrePol(n,:) - (n-1)*legendrePol(n-1,:))./n;
    end
end