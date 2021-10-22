function int = quasiMonteCarlo_ball(fun,N, dim, radius)
    arguments
        fun (:,:) double {mustBeNumeric, mustBeReal}
        N (1,1) double {mustBeNumeric, mustBeReal, mustBePositive}
        dim (1,1) double {mustBeNumeric, mustBeReal} = [] %mustBeFit2Fun(fun,dim,N)
        radius (1,1) double {mustBeNumeric, mustBeReal, mustBePositive} = 1
    end
    % Standard Quasi Monte Carlo Formulae
    if dim == 0
        int = sum(sum(fun))/N;
    else
        int = sum(fun,dim)/N;
    end
    % Standard Quasi Monte Carlo Formulae int = sum_i^N f(x_i) / N is for
    % inteegration over the s-dim cube [0,1]^s. The ball lies in [-radius,
    % radius]^s. Thus, volume element must be adapted. 
    int = int*(2*radius)^3;
end

% Verification of suitable dimensions costs a lot of time. Is dropped here
% for the sake of speedness
% function mustBeFit2Fun(f,dim,N)
%     if dim == 0
%         N_expected = length(f);
%     elseif dim == 1 || dim ==2 
%         N_expected = size(f,dim);
%     else
%         eid = 'Size:notEqual';
%         msg = 'dim is not expected. Must be {0,1,2}';
%         throwAsCaller(MException(eid,msg))
%     end
%     if N < N_expected
%         eid = 'Size:notEqual';
%         msg = 'Size of N is not expected';
%         throwAsCaller(MException(eid,msg))
%     elseif N > N_expected
%         warning('Number of summation points is smaller than N. If you are integrating over a subregion, its okay. Otherwise, check the code.');
%     end
% end