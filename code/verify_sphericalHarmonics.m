function  verify_sphericalHarmonics(N, phi, t, do_plot)
% verifies REAL fully normalized spherical harmonics

% N degree of spherical harmonics
% phi is the longitude; range phi in [0,2*pi[
% t is the polar distance; range t in [-1,1]
% can plot all spherical harmonics up to max degree N if do_plot is true

%function returns error if spherical harmonics are broken.
    arguments
        N  (1,1) double {mustBeNumeric,mustBeReal,mustBeInteger, mustBeSmaller(N)} 
        phi (:,1) double {mustBeNumeric,mustBeReal} 
        t (:,1) double {mustBeNumeric,mustBeReal,mustBeDimensionFit(phi,t)} 
        do_plot (1,1) logical = false 
    end
    
    % Generate evaluation grid in spherical and Cartesian coordinates.
    [Phi,T] = meshgrid(phi,t);
    x = sqrt(1-T.*T).*cos(Phi);
    y = sqrt(1-T.*T).*sin(Phi);
    z = T;
   
    %Note, that this function returnes real, fully normalized spherical
    %harmonics
    sphericalHarmonics_value = sphericalHarmonics(N, Phi, T);
    
    % Plots can be compared visually for instance [Michel 2013
    % 10.1007/978-0-8176-8403-7] p. 129
    if do_plot
        for i=1:size(sphericalHarmonics_value,1)
            figure;
            surf(x,y,z,squeeze(sphericalHarmonics_value(i,:,:))); 
            shading interp;
            colormap jet;
            [n,j] = transform_index2grad(i);
            title(['(', num2str(n), ',', num2str(j), ')'])
        end
    end
    
    % The real spherical harmoncs are calculated up to degree 4.
    sphericalHarmonicsTabular_value = zeros(size(sphericalHarmonics_value));
     if N >=0
        sphericalHarmonicsTabular_value(1,:,:) = 0.5/sqrt(pi);
    end
    if N >=1
        sphericalHarmonicsTabular_value(4,:,:) = sqrt(3/(4*pi)).*y;
        sphericalHarmonicsTabular_value(3,:,:) = sqrt(3/(4*pi)).*z;
        sphericalHarmonicsTabular_value(2,:,:) = sqrt(3/(4*pi)).*x;
    end
    if N >=2
        sphericalHarmonicsTabular_value(9,:,:) = 0.5*sqrt(15/pi).*x.*y;
        sphericalHarmonicsTabular_value(8,:,:) = 0.5*sqrt(15/pi).*y.*z;
        sphericalHarmonicsTabular_value(7,:,:) = 0.25*sqrt(5/pi).*(2*z.*z-x.*x-y.*y);
        sphericalHarmonicsTabular_value(6,:,:) = 0.5*sqrt(15/pi).*z.*x;
        sphericalHarmonicsTabular_value(5,:,:) = 0.25*sqrt(15/pi).*(x.*x-y.*y);
    end
    if N >=3
        sphericalHarmonicsTabular_value(16,:,:) = 0.25*sqrt(35/(2*pi)).*(3*x.*x-y.*y).*y;
        sphericalHarmonicsTabular_value(15,:,:) = 0.5*sqrt(105/pi).*x.*y.*z;
        sphericalHarmonicsTabular_value(14,:,:) = 0.25*sqrt(21/(2*pi)).*(4*z.*z-x.*x-y.*y).*y;
        sphericalHarmonicsTabular_value(13,:,:) = 0.25*sqrt(7/pi).*(2*z.*z-3*x.*x-3*y.*y).*z;
        sphericalHarmonicsTabular_value(12,:,:) = 0.25*sqrt(21/(2*pi)).*(4*z.*z-x.*x-y.*y).*x;
        sphericalHarmonicsTabular_value(11,:,:) = 0.25*sqrt(105/pi).*(x.*x-y.*y).*z;
        sphericalHarmonicsTabular_value(10,:,:) = 0.25*sqrt(35/(2*pi)).*(x.*x-3*y.*y).*x;
    end
    if N >= 4
        sphericalHarmonicsTabular_value(25,:,:) = 0.75*sqrt(35/pi).*(x.*x-y.*y).*x.*y;
        sphericalHarmonicsTabular_value(24,:,:) = 0.75*sqrt(35/(2*pi)).*(3*x.*x-y.*y).*z.*y;
        sphericalHarmonicsTabular_value(23,:,:) = 0.75*sqrt(5/pi).*(7*z.*z-1).*x.*y;
        sphericalHarmonicsTabular_value(22,:,:) = 0.75*sqrt(5/(2*pi)).*(7*z.*z-3).*z.*y;
        sphericalHarmonicsTabular_value(21,:,:) = (3/16)*sqrt(1/pi).*(35*z.^4-30*z.*z+3);
        sphericalHarmonicsTabular_value(20,:,:) = 0.75*sqrt(5/(2*pi)).*(7*z.*z-3).*x.*z;
        sphericalHarmonicsTabular_value(19,:,:) = (3/8)*sqrt(5/pi).*(x.*x-y.*y).*(7*z.*z-1);
        sphericalHarmonicsTabular_value(18,:,:) = 0.75*sqrt(35/(2*pi)).*(x.*x-3*y.*y).*x.*z;
        sphericalHarmonicsTabular_value(17,:,:) = (3/16)*sqrt(35/pi).*(x.*x.*(x.*x-3*y.*y) - y.*y.*(3*x.*x-y.*y));
    end
    
    %Calculate difference and through error if needed.
    err = max(max(max(abs(sphericalHarmonicsTabular_value - sphericalHarmonics_value))));
    if err > 1000*eps
        error('Spherical Harmonic Calculation is corruped.')
    end
end

function [n,j] = transform_index2grad(idx)
    n = ceil(sqrt(idx))-1;
    j = idx-n.^2-n-1;
end

function mustBeSmaller(a)
% Test for fitting dimension
    if a > 4
        eid = 'Size:notEqual';
        msg = 'N must be smaller than 4';
        throwAsCaller(MException(eid,msg))
    end
end

function mustBeDimensionFit(a,b)
% Test for fitting dimension
    if ~isequal(size(a),size(b)) 
        eid = 'Size:notEqual';
        msg = 'Size of t must be equal size of phi.';
        throwAsCaller(MException(eid,msg))
    end
end