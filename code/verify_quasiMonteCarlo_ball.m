function int = verify_quasiMonteCarlo_ball(testCase,N)
    arguments
        testCase (1,:) char {mustBeMember(testCase,{'one','x1_squared','x_squared_all','x_squared_vector'})} = 'x1_squared'
        N (1,1) double {mustBeNumeric, mustBeReal, mustBePositive} = 10000000
    end
    
    % Calculate Integration points via quasi-random sequence
    radius = 1;
    x = radius*calculate_integrationPoints(N);
    
    % Calculate real and approx integral for particular test case.
    switch testCase
        case 'one' 
            int_real = 4*pi*radius^3/3; % Unit ball integral 4*pi*/3
            fun = ones(size(x,1),1); % f(x) = 1
            int = quasiMonteCarlo_ball(fun,N,1,radius);
        case 'x1_squared'  
            int_real = 4*pi*radius^3/15; % Unit ball integral 4*pi*/15
            fun = x(:,1).^2; % f(x1,x2,x3) = x1^2
            int = quasiMonteCarlo_ball(fun,N,1,radius);
        case 'x_squared_all'
            int_real = 4*pi*radius^3/5; % Unit ball integral 4*pi*/5
            fun = [x(:,1).^2,x(:,2).^2,x(:,3).^2]; % f(x1,x2,x3) = x1^2
            int = quasiMonteCarlo_ball(fun,N,0,radius);
        case 'x_squared_vector'
            int_real = 4*pi*radius^3/15*ones(1,3); % Unit ball integral 4*pi*/5
            fun = [x(:,1).^2,x(:,2).^2,x(:,3).^2]; % f(x1,x2,x3) = x1^2
            int = quasiMonteCarlo_ball(fun,N,1,radius);
    end
    
    % Calculate difference, through error if needed.
    if max(max(abs((int_real-int)./int_real))) > 5e-05
       error('Quasi Monte Carlo Integration is corruped of to few points are chosen.')
    end
end