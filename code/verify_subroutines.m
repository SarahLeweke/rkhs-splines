function verify_subroutines(device,verify_closedSeriesRepresentations,verify_basisMathFunctions)
    arguments
        device (1,1) struct
        verify_closedSeriesRepresentations (1,1) logical = false
        verify_basisMathFunctions (1,1) logical = false
    end
    
    if verify_basisMathFunctions
        sequence = rand(1,7);
        verify_clenshawAlgorithm(6,repmat(sequence(1:7)',1,100),linspace(-1,1,100)) % For Legendre Polynomials with first and second derivative up to degree 6
        verify_sphericalHarmonics(4,linspace(0,2*pi,100),linspace(-1,1,100)); % Spherical Harmonics up to degree 4
        verify_quasiMonteCarlo_ball
    end
    
    if verify_closedSeriesRepresentations
        %Verification of closed integral kernel nabla_y K_M(x,y) 
        verify_closedRepresentation_IntegralKernel(rand(1000,3));
        %Verify Spline Kernel
        x = rand(1000,3); x = device.radiusCerebrum.*x./max(sqrt(sum(x.^2,2)));
        z = rand(1,3); z = device.radiusCerebrum.*z./sqrt(sum(z.^2,2));
        h = 0.8;
        verify_closedRepresentation_SplineKernel(x,z,1000,device.radiusCerebrum,h);
        verify_closedRepresentation_SplineKernel_laplace(x,z,1000,device.radiusCerebrum,h);
        %Verification of quasi Monte Carlo integration over a ball. 
    end 
end