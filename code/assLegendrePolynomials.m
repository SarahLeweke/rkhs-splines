function assLeg = assLegendrePolynomials(N, t)
% Leweke 2021
% This algorithm computes the associated legendre functions with the
% constant c_{n,k} = sqrt((2n+1)/(4pi) * (n-k)!/(n+k)!)

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

    % Initialize   
    assLeg = zeros(N+1,N+1, length(t));
    assLeg(1,1,:) = ones(length(t),1)./sqrt(4*pi);

    % algorithm can be found in for instance [Fengler Diss 2005 U Kaiserslautern Algo 15.3.3]
    for n=1:N
        assLeg(n+1,n+1,:) = sqrt((2*n+1)/(2*n))*sqrt(1-t.^2).*squeeze(assLeg(n,n,:))';
    end
    for k=0:N-1
        assLeg(k+2, k+1,:) = sqrt(2*k+3)*t.*squeeze(assLeg(k+1, k+1, :))';
        for n=(k+2):N
            assLeg(n+1,k+1,:) = sqrt(((2*n-1)*(2*n+1))/((n-k)*(n+k)))*t.*squeeze(assLeg(n,k+1,:))' ...
                - sqrt((2*n+1)/(2*n-3))*sqrt(((n+k-1)*(n-k-1))/((n-k)*(n+k)))*squeeze(assLeg(n-1,k+1,:))';
        end
    end
end