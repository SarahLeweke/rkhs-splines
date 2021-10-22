function splineMat_closed = integragte_quasiMonteCarlo_ball_scalarMEG(device,h,Nx,Nz, num_diff,parallel_computing)
    arguments
        device (1,1) struct
        h (1,1) double {mustBeNumeric, mustBeReal, mustBePositive} = 0.8
        Nx (1,1) double {mustBeNumeric, mustBeReal, mustBePositive, mustBeInteger} = 4500000
        Nz (1,1) double {mustBeNumeric, mustBeReal, mustBePositive} = 1*Nx
        num_diff (1,1) logical = true
        parallel_computing (1,1) logical = false
    end
    
    x = device.radiusCerebrum*calculate_integrationPoints(Nx);
    kernel_xy = calculate_ingetralKernel_ScalarMEG(x,device.MEG.measurementPosition,device.MEG.normals);
    if Nx == Nz
        z = x;
        kernel_zy = kernel_xy;
    else
        z = device.radiusCerebrum*calculate_integrationPoints(Nz);
        kernel_zy = calculate_ingetralKernel_ScalarMEG(z,device.MEG.measurementPosition,device.MEG.normals);
    end
    
    % Calculates a closed representation of the scalar MEG integral kernel
    
    % Time and memory resources are split. One dimension (z) is
    % integrated via for loop, the other (x) via vectorization
    splineMat_closed = zeros(device.MEG.numberMeasurement);

    if parallel_computing
        parfor i=1:size(z,1)
            spline_term = calculate_functional2SplineKernel_ScalarMEG_laplace(x,z(i,:),device.radiusCerebrum, h, num_diff,8)';
            tmp = (2*device.radiusCerebrum)^3*spline_term*kernel_xy/Nx; %quasiMonteCarlo_ball(spline_term*kernel_xy,Nx, 1, device.radiusCerebrum)';
             % Vectorized variant of quasi-Monte-Carlo integration, Scaling
             % factor for ball integration is multiplied later
            splineMat_closed = splineMat_closed + (tmp'*kernel_zy(i,:))/Nz;
            if mod(i,1000) == 0
               disp([num2str(100*i/size(z,1)), '% calculated.']) 
            end
        end
    else
        for i=1:size(z,1)
            spline_term = calculate_functional2SplineKernel_ScalarMEG_laplace(x,z(i,:),device.radiusCerebrum, h, num_diff,false)';
            tmp = (2*device.radiusCerebrum)^3*spline_term*kernel_xy/Nx; %quasiMonteCarlo_ball(spline_term*kernel_xy,Nx, 1, device.radiusCerebrum)';
             % Vectorized variant of quasi-Monte-Carlo integration, Scaling
             % factor for ball integration is multiplied later
            splineMat_closed = splineMat_closed + (tmp'*kernel_zy(i,:))/Nz;
        end
    end
    % Scaling factor for quasi Monte-Carlo on balls and transformation to fT
    splineMat_closed = splineMat_closed*(2*device.radiusCerebrum)^3*1e30*device.mu0^2;
end
