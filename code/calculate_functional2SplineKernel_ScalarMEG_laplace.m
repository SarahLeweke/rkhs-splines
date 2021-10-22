function fun = calculate_functional2SplineKernel_ScalarMEG_laplace(x,z,rho, kernel_h, num_diff,diff_order,h, threshold)
    arguments
        x (:,3) double {mustBeNumeric, mustBeReal}
        z (:,3) double  {mustBeNumeric, mustBeReal}
        rho (1,1) double {mustBeNumeric, mustBeReal} = 0.071
        kernel_h (1,1) double {mustBeNumeric, mustBeReal} =  0.8
        num_diff (1,1) logical = true
        diff_order (1,1) double {mustBeNumeric, mustBePositive, mustBeInteger} = 8
        h (1,1) double {mustBeNumeric, mustBeReal} = 3e-3; %7.5e-4;
        threshold (1,1) double {mustBeNumeric, mustBeReal} = 1e-2
    end
    
    if num_diff 
    % Calculates the three-dimensional Laplacian for the function obtained
    % by calculate_functional2SplineKernel_ScalarMEG with respect to the
    % desired accuracy order
   
        f = calculate_functional2SplineKernel_ScalarMEG( x, z, rho, kernel_h);
  
        fun = zeros(size(f));
        
        if diff_order == 2
            for i = 1:size(x, 2)
                ind = (abs(x(:,i)) >= threshold);
                hr = h .* (abs(x(:, i)) .* ind + (~ind));
                m = zeros(size(x));
                m(:,i) = hr;

                left = calculate_functional2SplineKernel_ScalarMEG( x - m, z, rho, kernel_h);
                right = calculate_functional2SplineKernel_ScalarMEG( x + m, z, rho, kernel_h);

                fun = fun + (right - 2 .* f + left) ./ repmat((hr.^2).', size(z, 1), 1);
            end
        elseif diff_order == 6
            for i = 1:size(x, 2)
                ind = (abs(x(:,i)) >= threshold);
                hr = h .* (abs(x(:, i)) .* ind + (~ind));
                m = zeros(size(x));
                m(:,i) = hr;

                left = calculate_functional2SplineKernel_ScalarMEG( x - m, z, rho, kernel_h);
                left2 = calculate_functional2SplineKernel_ScalarMEG( x - 2*m, z, rho, kernel_h);
                left3 = calculate_functional2SplineKernel_ScalarMEG( x - 3*m, z, rho, kernel_h);
                right = calculate_functional2SplineKernel_ScalarMEG( x + m, z, rho, kernel_h);
                right2 = calculate_functional2SplineKernel_ScalarMEG( x + 2*m, z, rho, kernel_h);
                right3 = calculate_functional2SplineKernel_ScalarMEG( x + 3*m, z, rho, kernel_h);

                fun = fun + (2*left3-27*left2+270*left-490*f+270*right-27*right2+2*right3) ./ repmat(180*(hr.^2).', size(z, 1), 1);
            end
        elseif diff_order == 8
            for i = 1:size(x, 2)
                ind = (abs(x(:,i)) >= threshold);
                hr = h .* (abs(x(:, i)) .* ind + (~ind));
                m = zeros(size(x));
                m(:,i) = hr;

                left = calculate_functional2SplineKernel_ScalarMEG( x - m, z, rho, kernel_h);
                left2 = calculate_functional2SplineKernel_ScalarMEG( x - 2*m, z, rho, kernel_h);
                left3 = calculate_functional2SplineKernel_ScalarMEG( x - 3*m, z, rho, kernel_h);
                left4 = calculate_functional2SplineKernel_ScalarMEG( x - 4*m, z, rho, kernel_h);
                right = calculate_functional2SplineKernel_ScalarMEG( x + m, z, rho, kernel_h);
                right2 = calculate_functional2SplineKernel_ScalarMEG( x + 2*m, z, rho, kernel_h);
                right3 = calculate_functional2SplineKernel_ScalarMEG( x + 3*m, z, rho, kernel_h);
                right4 = calculate_functional2SplineKernel_ScalarMEG( x + 4*m, z, rho, kernel_h);

                fun = fun + (-left4/560 + 8*left3/315 - left2/5 + 8*left/5 - 205*f/72 + 8*right/5 - right2/5 + 8*right3/315 - right4/560) ./ repmat((hr.^2).', size(z, 1), 1);
            end
       
        end
    fun = fun';
        
    else
        % Calculates the Laplacian by means of the analytic derivative of
        % the Legendre series and the Clenshaw summation. Method is more
        % inaccuarate due to summation and significant slower than
        % numerical differentation in our test cases.
        r = sqrt(sum(x.^2,2));
        xi = x./r;
        fun = zeros(size(x,1),size(z,1));
        N = 1000;
        n = linspace(0,N,N+1);
           
        for i=1:size(z,1)
            eta = z(i,:)./norm(z(i,:),2);

            % Generate summation coefficients.
            coeff = (kernel_h*kernel_h*r.*norm(z(i,:),2)/rho^2).^n.*(2*n+3).^2;

            % Calculate series representation via Clenshaw Algorithm
            fun(:,i) = kernel_h*kernel_h.*clenshawLegendre(coeff',xi*eta')'/(pi*rho^5);
        end
    end
    
    % Besides these two methods we have also tested automatic
    % differentiation from the ADIMAT package. Since it could not cover the
    % vectorized structure of the code, it was too slow to use it for the
    % desire amount of evaluation points.
%     adOpts = admOptions();
%     adOpts.hessianStrategy = 'for2';
%     adOpts.independents = 1;
%     H = admHessian(@calculate_functional2SplineKernel_ScalarMEG_2, 1, x, z, rho, h, adOpts);
% 
%     nx = size(x, 1);
%     nz = size(z, 1);
%     fun = zeros(nz, nx);
%     k = 1;
%     for j = 1:size(fun,2)
%         for i = 1:size(fun,1)
%             hh = H([1:nx:3*nx]+j-1, [1:nx:3*nx]+j-1,k);
%             fun(i,j) = trace(hh);
%             k = k + 1;
%         end
%     end
    
end