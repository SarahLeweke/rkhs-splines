function verify_closedRepresentation_SplineKernel(x,z,N,rho,h)
    arguments
        x (:,3) double {mustBeNumeric, mustBeReal}
        z (1,3) double {mustBeNumeric, mustBeReal} = 10*x(4,:)/norm(x(4,:),2)
        N (1,1) double {mustBeNumeric, mustBeReal, mustBePositive} = 1000
        rho (1,1) double {mustBeNumeric, mustBeReal} = 0.071
        h (1,1) double {mustBeNumeric, mustBeReal} =  0.9
    end
    
    % Calculate radius and angular. 
    r = sqrt(sum(x.^2,2));
    xi = x./r;
    eta = z./norm(z,2);
    
    % Generate summation coefficients.
    n = linspace(0,N,N+1);
    coeff = (h*h*r.*norm(z,2)/rho^2).^n.*(2*n+3);
    %coeff(:,1) = zeros(1,size(xi,1));
    
    % Calculate series representation via Clenshaw Algorithm
    integralKernel_series = (h*r).^2.*clenshawLegendre(coeff',xi*eta')'/(2*pi*rho^5);
   
    % Calculate closed representation 
    integralKernel_close = calculate_functional2SplineKernel_ScalarMEG(x,z,rho, h)';
 
   
    % Calculate difference and through error if needed.
    if max(max(abs((integralKernel_close-integralKernel_series)./integralKernel_series))) > 8e-15
        disp('Integral kernel represenation is broken.')
    end 
end