function fun = calculate_ingetralKernel_ScalarMEG(x,y,nu) 
    arguments
        x (:,3) double {mustBeNumeric, mustBeReal}
        y (:,3) double {mustBeNumeric, mustBeReal} = 10*x(1,:)/norm(x(1,:),2)
        nu (:,3) double {mustBeNumeric, mustBeReal,mustBeFit(y,nu)} = y
    end
    
    % Formula is stated in Sec. 5.1 of [3] in the itemize over avoide
    % inverse crime
    
    r = sqrt(sum(x.^2,2));
    normDiff = zeros(size(x,1),size(y,1));
    for k=1:size(y,1)
        normDiff(:,k) = sqrt(sum((x-y(k,:)).^2,2));
    end
    s = sqrt(sum(y.^2,2));
    xi = x./r;
    eta = y./s;
    t = xi*eta';
   
    idx = abs(t - 1) <1e-13;
    idx2 = abs(t + 1) <1e-13;
    
    fun =  zeros(size(x,1),size(y,1));
    
    % Due to the singularity in the case of parallel directions, three
    % cases need to be considered. 
    if sum(sum(idx)) >= 1
        tmp = (1./(ones(size(x,1),1)*s').^2 - 1./((ones(size(x,1),1)*s').*(ones(size(x,1),1)*s'-r*ones(1,size(y,1))))).*(ones(size(x,1),1)*sum(nu.*eta,2)');
        fun(idx) = tmp(idx);
    end
    if sum(sum(idx2)) >= 1
        tmp =  -(-1./(ones(size(x,1),1)*s').^2 + 1./((ones(size(x,1),1)*s').*((ones(size(x,1),1)*s')+r*ones(1,size(y,1))))).*(ones(size(x,1),1)*sum(nu.*eta,2)');
        fun(idx2) = tmp(idx2);
    end
    
    idxRest = ~(idx+idx2);
    if sum(idxRest) >= 1
        tmp = (((ones(size(x,1),1)*s')-r*ones(1,size(y,1)).*t)./normDiff-1).*(xi*nu') - (((ones(size(x,1),1)*s').*t-r*ones(1,size(y,1)))./normDiff - t).*(ones(size(x,1),1)*sum(nu.*eta,2)');
        tmp = tmp./(r*ones(1,size(y,1)).*(t.^2-1).*(ones(size(x,1),1)*s')) + (ones(size(x,1),1)*sum(nu.*eta,2)')./((ones(size(x,1),1)*s').*(ones(size(x,1),1)*s'));
        
        fun(idxRest) = tmp(idxRest);
    end
    fun = fun./(4*pi);
end

function mustBeFit(a,b)
    if ~isequal(size(a,1),size(b,1))
        eid = 'Size:notEqual';
        msg = 'y and nu must have the same amount of points';
        throwAsCaller(MException(eid,msg))
    end
end