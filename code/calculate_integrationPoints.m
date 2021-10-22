function points = calculate_integrationPoints(N,seed,g)
    arguments
       N (1,1) double {mustBeNumeric, mustBeReal, mustBePositive} = 2000;
       seed (1,1) double {mustBeNumeric, mustBeReal, mustBePositive} = 0.5
       g (1,1) double {mustBeNumeric, mustBeReal}  = 1.22074408460575947536
    end
    % Calculates Quasi-Monte Carlo integration points based on quasi-random
    % sequence.
    % More information can be found in 
    % http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/
    
    a1 = 1.0/g;
    a2 = 1.0/(g*g);
    a3 = 1.0/(g*g*g);

    n = 1:N;
    x = mod(seed+a1*n,1); 
    y = mod(seed+a2*n,1);
    z = mod(seed+a3*n,1); 
    
    points = 2*([x',y',z']-0.5);
    r = sqrt(sum(points.^2,2));
    points = points(r<=1,:);
end