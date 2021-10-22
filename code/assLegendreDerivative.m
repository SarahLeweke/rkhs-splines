function [assLegendreDeriv, assLegendre] = assLegendreDerivative(N, t)
% Leweke 2021
% This algorithm computes the first derivatives of the function from
% assLegendrePolynomials, i.e. from the associated legendre Polynomials
% with the normalization constant  c_{n,k} = sqrt((2n+1)/(4pi) * (n-k)!/(n+k)!)

% N maximal degree of corresponding Legendre polynomials
% t vector of evaluation points
% range t in [-1,1]

% returns assLeg with dim(assLeg) = (N+1) x (N+1) x length(t)
% first dimension = degree of Legendre polynomials
% second dimension = order of ass. Legendre functions
    arguments
        N  (1,1) double {mustBeNumeric,mustBeReal,mustBeInteger} 
        t (:,1) double {mustBeNumeric,mustBeReal} 
    end

    % Linearize t (row vector)
    t = t';
    
    % Compute the polynomials
    assLegendre = assLegendrePolynomials(N, t); % size : (N+1)x(N+1)x(length(cosPhi))
    assLegendreDeriv = zeros(N+1, N+1, length(t));
    
    legendreD = legendreDeriv(N,t);
    
    % for k=0
    Nhelp = linspace(0,N,N+1)';
    Nhelp = repmat(Nhelp,1,length(t));
    assLegendreDeriv(:,1,:) = legendreD.*sqrt((2*Nhelp+1)./(4*pi));
    
    % algorithm can be found in for instance [Fengler Diss 2005 U Kaiserslautern Algo 15.3.6]
    for n= 1:N
       % for k=n 
       idx = find(abs(abs(t)-1) > 1e-5);
       if ~isempty(idx)
            assLegendreDeriv(n+1,n+1,idx) = -n*t(idx)./(1-t(idx).^2).*squeeze(assLegendre(n+1,n+1,idx))';
       end
       idx2 = find(abs(abs(t)-1) <= 1e-5);
       if n == 0 || n > 2
           assLegendreDeriv(n+1,n+1,idx2) = zeros(1,length(t(idx2)));
       elseif n == 1
           assLegendreDeriv(n+1,n+1,idx2) = Inf;
       elseif n == 2
           assLegendreDeriv(n+1,n+1,idx2) = (-1)*sqrt((2*n+1)/(4*pi*4*3*2))*6.*sign(t(idx2));
       end
       % for the other k's    
       for k = 1:(n-1)
           if ~isempty(idx)
                assLegendreDeriv(n+1,k+1,idx) = (-k*t(idx))./(1-t(idx).^2).*squeeze(assLegendre(n+1,k+1,idx))' + ...
                sqrt(((n+k+1)*(n-k))./(1-t(idx).^2)).*squeeze(assLegendre(n+1,k+2,idx))';
           end
           if ~isempty(idx2)
                if k == 1 
                    assLegendreDeriv(n+1,k+1,idx2) =  t(idx2).*legendreD(n+1,idx2)./sqrt(1-t(idx2).^2);
                elseif k == 2 || n >= 2
                    legendre2D = 1;
                    for j=0:floor(n/2)
                        legendre2D = legendre2D*(1 + (-1)*(n-2*j-3)*(n-2*j-2)/(2*(2*n-2*j-1)*(j+1)));
                    end
                    assLegendreDeriv(n+1,k+1,idx2) = -2*t(idx2).*sign(t(idx2)).^(n).*legendre2D; %This is the correct solution
                else 
                    assLegendreDeriv(n+1,k+1,idx2) = zeros(1,length(t(idx2)));
                end
           end
       end
    end
end