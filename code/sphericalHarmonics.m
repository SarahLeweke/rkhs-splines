function spH = sphericalHarmonics(N, phi, t)
% Leweke 2021

% computes REAL fully normalized spherical harmonics

% N degree of spherical harmonics
% phi is the longitude; range phi in [0,2*pi[
% t is the polar distance; range t in [-1,1]

% returns spH of size (N+1)^2 x length(t)
% spH are grouped row-wise of order (0,0), (1,-1), (1,0), (1,-1), (2,-2)...
    arguments
        N  (1,1) double {mustBeNumeric,mustBeReal,mustBeInteger} 
        phi (:,:) double {mustBeNumeric,mustBeReal} 
        t (:,:) double {mustBeNumeric,mustBeReal,mustBeDimensionFit(phi,t)} 
    end
    
    dim = size(phi);
    % Linearize vectors
    t = t(:);
    phi = phi(:);

    leg = assLegendrePolynomials(N, t');
    spH = zeros((N+1)^2, length(t));
    
    % formula for real spherical harmonics for instance [Michel 2013
    % 10.1007/978-0-8176-8403-7] Eq. (5.161)
    spH(1,:) = squeeze(leg(1,1,:));
    for n = 1:N
       spH(n^2+n+1,:) = squeeze(leg(n+1,1,:));
       j=(1:n)';
       counter = n^2+n+1;
       spH(counter-j,:) = shiftdim(leg(n+1,j+1,:),1).*cos(j.*phi')*sqrt(2);
       spH(counter+j,:) = shiftdim(leg(n+1,j+1,:),1).*sin(j.*phi')*sqrt(2);
       % No need of changes for j = 0
    end
    
    if min(dim)>1
        spH = reshape(spH,size(spH,1),dim(1),dim(2));
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