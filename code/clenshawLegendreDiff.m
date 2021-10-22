function sum = clenshawLegendreDiff(coeff, t)
% Leweke 2021
% Calculates sum of Legendre polynomials up to degree length(coeff,1). 
% size(coeff) = degree x evaluation points || degree x 1 for spatially
% independent coefficients
% size(t) = evaluation points x 1
% range(t) in [-1,1]

% Code evaluates the summation at several points
    arguments
        coeff (:,:) double {mustBeNumeric,mustBeReal} 
        t (:,1) double {mustBeNumeric,mustBeReal,mustBeDimensionFit(coeff,t)} 
    end
    
    % Initialization
    t = t(:);
    N = size(coeff,1)-1;   
    U1 = zeros(length(t),1);
    U2 = zeros(length(t),1);
 
    
    for k=(N-1):-1:1
       % algorithm can be found in for instance [Fengler Diss 2005 U Kaiserslautern Algo 15.2.2]
       U0 = (2*k+3)/(k+1)*t.*U1 - (k+3)/(k+2)*U2 + coeff(k+2,:)';
       U2 = U1;
       U1 = U0;
    end
    sum = coeff(2,:) - 1.5*U2' + (3*U1.*t)';
end

function mustBeDimensionFit(a,b)
% Test for fitting dimension
    if ~isequal(size(a,2),length(b(:))) && ~isequal(size(a,2),1)
        eid = 'Size:notEqual';
        msg = 'Second dimension of first input must equal length of second input.';
        throwAsCaller(MException(eid,msg))
    end
end