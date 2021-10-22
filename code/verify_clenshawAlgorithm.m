function verify_clenshawAlgorithm(N, coeff, t)
% Leweke 2021
% Checks up to max degree 6 if Clenshaw's Algorithms returns correct
% results.
% size(coeff) = degree x evaluation points || degree x 1 for spatially
% independent coefficients
% size(t) = evaluation points x 1
% range(t) in [-1,1]

%function returns error if Clenshaw is broken.
    arguments
        N  (1,1) double {mustBeNumeric,mustBeReal, mustBeSmaller(N)} 
        coeff (:,:) double {mustBeNumeric,mustBeReal, mustBeDimensionFitN(coeff,N)} 
        t (:,1) double {mustBeNumeric,mustBeReal,mustBeDimensionFit(coeff,t)} 
    end
    
    % Implement closed representation of Legendre polynomials and its first
    % two derivatives. See https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas
    p = zeros(N+1,length(t));
    pDiff = zeros(N+1,length(t));
    pDiff2 = zeros(N+1,length(t));
    if N >= 0
            p(1,:) = ones(size(t))';
    end
    if N >= 1
            p(2,:) = t';
            pDiff(2,:) = ones(size(t))';
    end
    if N >= 2
            p(3,:) = (0.5*(3*t.*t -1))';
            pDiff(3,:) = 3*t';
            pDiff2(3,:) = 3*ones(size(t))';
    end
    if N >= 3
            p(4,:) = (0.5*(5*t.*t.*t-3*t))';
            pDiff(4,:) = 0.5*(15*t.*t - 3);
            pDiff2(4,:) = 15*t';
    end
    if N >= 4
            p(5,:) = 1/8*(35*t.^4 - 30*t.*t+3)';
            pDiff(5,:) = 1/8*(4*35*t.^3-60*t)';
            pDiff2(5,:) = 1/8*(12*35*t.^2-60)';
    end
    if N >= 5
            p(6,:) = (1/8*(63*t.^5 - 70*t.^3 + 15*t))';
            pDiff(6,:) = 1/8*(5*63*t.^4-210*t.^2+15)';
            pDiff2(6,:) = 1/8*(20*63*t.^3-420*t)';
    end
    if N >= 6
            p(7,:) = (1/16*(231*t.^6-315*t.^4+105*t.*t-5))';
            pDiff(7,:) = 1/16*(6*231*t.^5-4*315*t.^3+210*t)';
            pDiff2(7,:) = 1/16*(30*231*t.^4-12*315*t.^2+210)';
    end
    
    % Use standard summation for explicit representation. Summation errors
    % are neglected since maximal degree is 6.
    sumLeg = sum(coeff.*p,1);
    sumDiff = sum(coeff.*pDiff,1);
    sumDiff2 = sum(coeff.*pDiff2,1);
   
    % Calculate summation via Clenshaw Algorithm
    clenshaw_sum = clenshawLegendre(coeff,t);
    clenshaw_sumDiff = clenshawLegendreDiff(coeff,t);
    clenshaw_sumDiff2 = clenshawLegendreDiff2(coeff,t);
    
    % Calculate difference and through error if needed.
    err = max(abs(sumLeg - clenshaw_sum)./sumLeg);
    errDiff = max(abs(sumDiff - clenshaw_sumDiff)./sumDiff);
    errDiff2 = max(abs(sumDiff2 - clenshaw_sumDiff2)./sumDiff2);
    if err > 1000*eps
        disp(num2str(err))
        error('Clenshaw Algorithm is corruped.')
    end
    if errDiff > 1000*eps
       error('Clenshaw Algorithm for first Derivative is corruped.')
    end
    if errDiff2 > 2000*eps
       error('Clenshaw Algorithm for second Derivative is corruped.')
    end
end

function mustBeSmaller(a)
% Make sure all closed representation up to maximum degree are implemented.
    if a > 6
        eid = 'Size:notEqual';
        msg = 'N must be smaller than 7';
        throwAsCaller(MException(eid,msg))
    end
end

function mustBeDimensionFitN(a,b)
% Test for fitting dimension
  if ~isequal(size(a,1),b+1) 
        eid = 'Size:notEqual';
        msg = 'First dimension of coeff  must equal N+1.';
        throwAsCaller(MException(eid,msg))
  end
end

function mustBeDimensionFit(a,b)
% Test for fitting dimension
    if ~isequal(size(a,2),length(b(:))) && ~isequal(size(a,2),1)
        eid = 'Size:notEqual';
        msg = 'Second dimension of first input must equal length of second input.';
        throwAsCaller(MException(eid,msg))
    end
end