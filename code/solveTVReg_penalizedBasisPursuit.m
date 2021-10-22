function x = solveTVReg_penalizedBasisPursuit(g,M,lambda, x0, epsilon, tau, sigma,  z)
    % Algorithm based on Functional penalised basis pursuit on spheres [Simeoni 2021]
    % https://doi.org/10.1016/j.acha.2020.12.004
    arguments
        g (:,1)  {mustBeNumeric, mustBeReal} 
        M (:,:)  {mustBeNumeric, mustBeReal} 
        lambda (1,1)  {mustBeNumeric, mustBeReal, mustBePositive} 
        x0 (:,1) {mustBeNumeric, mustBeReal, mustBeSizeFit(g,M,x0)} = zeros(size(g))
        epsilon (1,1)  {mustBeNumeric, mustBeReal, mustBePositive} = 100*eps
        tau (1,1)  {mustBeNumeric, mustBeReal, mustBePositive, mustBePositive} = 1/norm(M)
        sigma   (1,1)  {mustBeNumeric, mustBeReal, mustBePositive} = 1/norm(M)
        z (:,1) {mustBeNumeric, mustBeReal, mustBeSameSize(x0,z)} = zeros(size(g))
    end
    delta = 0;
    n = 0;
    
    while delta == 0 && n <10000
        n = n + 1;
        x  = soft(lambda, tau, x0 - tau*(M*z));
        z = proxCon(sigma,g,z + sigma*M*(2*x-x0));
        if norm(x-x0,2) < epsilon*norm(x0,2)
            delta = 1;
        else
            x0 = x;
        end
    end
 end

function x = soft(lambda, tau, x)
    arguments
        lambda (1,1)  {mustBeNumeric, mustBeReal, mustBePositive} 
        tau (1,1)  {mustBeNumeric, mustBeReal, mustBePositive} 
        x (:,1) {mustBeNumeric, mustBeReal}
    end

    x = max(abs(x)-tau*lambda,0).*sign(x);
end

function z = proxCon(sigma,g,z)
    arguments
        sigma   (1,1)  {mustBeNumeric, mustBeReal, mustBePositive} 
        g (:,1)  {mustBeNumeric, mustBeReal} 
        z (:,1) {mustBeNumeric, mustBeReal, mustBeSameSize(z,g)}
    end

    z = z - sigma*prox_L2ball(g,z/sigma); %tau = 1/sigma, aber prox_L2ball unabhängig von tau
end

function z = prox_L2ball(g,z,noiseLevel)
    arguments
        g (:,1)    {mustBeNumeric, mustBeReal} 
        z (:,1) {mustBeNumeric, mustBeReal, mustBeSameSize(z,g)}
        noiseLevel (1,1)   {mustBeNumeric, mustBeReal} = 1e-6
    end
    
    %proximal operator of i_(B2,noiseLevel)(Mx-g) mit 
    z = noiseLevel*(z - g)/(norm(z-g,2)) + g;
    
end

function mustBeSizeFit(a,b,c)
% Test for fitting dimension
    if ~isequal(length(a),length(b)) || ~isequal(size(b,1),size(b,2)) ||~isequal(length(a),length(c)) 
        eid = 'Size:notEqual';
        msg = 'Data, splineMatrix and coefficient vector must be dimensional fit dim(g) = dim (M x)';
        throwAsCaller(MException(eid,msg))
    end
end

function mustBeSameSize(a,b)
    if ~isequal(length(a), length(b)) 
        eid = 'Size:notEqual';
        msg = 'x and z must have same size';
        throwAsCaller(MException(eid,msg))
    end
end