function fun = calculate_functional2SplineKernel_ScalarMEG(x,z,rho, h)
    arguments
        x (:,3) double {mustBeNumeric, mustBeReal}
        z (:,3) double  {mustBeNumeric, mustBeReal}
        rho (1,1) double {mustBeNumeric, mustBeReal} = 0.071
        h (1,1) double {mustBeNumeric, mustBeReal} =  0.9
    end
    % calculates |x| Laplace_z |z| K(x,z)  with K(x,z) reproducing kernel from
    % spline
    % Formulae is stated in [3], sec. 5.1., itemize over avoid inverse
    % crime
    t = rho/h;
    z_rad = sqrt(sum(z.^2,2));
    x_scaled = x./t;
    x_scaled_rad = sqrt(sum(x_scaled.^2,2));
    fun = (3*t^2-4*z*x' + z_rad.^2*x_scaled_rad'.^2)./((z_rad.^3*ones(size(x_scaled_rad'))).*(ones(size(z_rad))*x_scaled_rad'.^2+(t*t-2*z*x')./(z_rad.^2*ones(size(x_scaled_rad')))).^(3/2));  
    fun = fun.*(ones(size(z_rad))*(x_scaled_rad.*sqrt(sum(x.^2,2)))')/(2*pi*rho^3);
end