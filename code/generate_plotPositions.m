function [eta_x,eta_y,eta_z,plotRadius,Phi,T] = generate_plotPositions(device, N, plotRadiusPerCent)
    arguments
        device (1,1) struct
        N (1,1) double {mustBeNumeric, mustBeReal} = 60
        plotRadiusPerCent (1,1) double {mustBeNumeric, mustBeReal} = 0.99
    end
    plotRadius = device.radiusCerebrum*plotRadiusPerCent;
        
    % Check if plot positions are already calculated - for speedup
    if exist('../solution/plot_positions.mat','file') == 2
        plotRadiusPerCent_recent = plotRadiusPerCent;
        load('../solution/plot_positions.mat', 'Phi', 'plotRadiusPerCent','x','y','z','T');
        if length(Phi) ~= N || plotRadiusPerCent_recent ~= plotRadiusPerCent
            delete('../Solution/plot_positions.mat');
        else
            eta_x = x./plotRadius;
            eta_y = y./plotRadius;
            eta_z = z./plotRadius;
        end
    end 
    
    % If necessary, calculate plot positions (euklidean and spherical coordinates) as points on an equiangular Driscoll-Healy
    % grid
    if ~exist('../solution/plot_positions.mat','file')
        t = linspace(-1,1,N);
        phi = linspace(0,2*pi,N);

        [Phi,T] = meshgrid(phi,t);

        eta_x = sqrt(1-T.*T).*cos(Phi);
        eta_y = sqrt(1-T.*T).*sin(Phi);
        eta_z = T;
        
        plotFile = struct('type', 'scalarMEG',  'x', plotRadius*eta_x, 'y', plotRadius*eta_y, 'z', plotRadius*eta_z,'plotRadiusPerCent', plotRadiusPerCent, 'Phi', Phi, 'T',T);
        save('../solution/plot_positions.mat', '-struct', 'plotFile'); 
    end
end