function verify_closedRepresentation_SplineKernel_laplace(x,z,N,rho,h)
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
    coeff = (h*h*r.*norm(z,2)/rho^2).^n.*(2*n+3).^2;

    
    % Calculate series representation via Clenshaw Algorithm
    integralKernel_series = h*h.*clenshawLegendre(coeff',xi*eta')'/(pi*rho^5);
   
    % Calculate closed representation 
    integralKernel_close = calculate_functional2SplineKernel_ScalarMEG_laplace(x,z,rho, h,true, 8, 1e-3,1e-5)';
 
   
    % Calculate difference and through error if needed.
    if max(max(abs((integralKernel_close-integralKernel_series)./integralKernel_series))) > 7.5e-5
        disp('Laplace of Integral kernel represenation is broken.')
        disp(max(max(abs((integralKernel_close-integralKernel_series)./integralKernel_series))))
    end
end