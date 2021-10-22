function verify_closedRepresentation_IntegralKernel(x,y,N)
    arguments
        x (:,3) double {mustBeNumeric, mustBeReal}
        y (1,3) double {mustBeNumeric, mustBeReal} = 10*x(4,:)/norm(x(4,:),2)
        N (1,1) double {mustBeNumeric, mustBeReal, mustBePositive} = 100
    end
    
    % Calculate radius and angular. 
    r = sqrt(sum(x.^2,2));
    xi = x./r;
    eta = y./norm(y,2);
    
    % Generate summation coefficients.
    k = linspace(0,N,N+1);
    coeff = (r./norm(y,2)).^k;
    coeff(:,1) = zeros(1,size(xi,1));
    coeffDiff = coeff./(k+1);
    coeffDiff(:,1) = zeros(1,size(xi,1));
    
    % Calculate series representation via Clenshaw Algorithm
    integralKernel_series = -eta.*clenshawLegendre(coeff',xi*eta')' ...
        + (xi - (xi*eta').*eta).*clenshawLegendreDiff(coeffDiff',xi*eta')';
    integralKernel_series = integralKernel_series/(4*pi*norm(y,2)^2);
   
    % Calculate closed representation 
    integralKernel_close =  calculate_ingetralKernel_ScalarMEG(x,y);
   
    % Calculate difference and through error if needed.
    if max(max(abs((integralKernel_close-sum(integralKernel_series*y',2))./sum(integralKernel_series*y',2)))) > 5e-10
        error('Integral kernel represenation is broken.')
    end
end