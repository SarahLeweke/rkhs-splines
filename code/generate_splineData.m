function [data, caseParameter] = generate_splineData(device, data_generatingCase, prob_case, plot_program, visualize_current)
    arguments
        device (1,1) struct
        data_generatingCase (1,:) char {mustBeMember(data_generatingCase,{'const','onb','2splines','2splines_avoidIC', '2splines_vec'})} = 'const'
        prob_case (1,:) char {mustBeMember(prob_case, {'scalarMEG', 'vectorMEG', 'vectorEEG'})} = 'scalarMEG'
        plot_program (1,:) char {mustBeMember(plot_program, {'none','matlab', 'tikz'})} = 'tikz'
        visualize_current (1,1) logical = false
    end
    
    % Set all parameters corresponding to the chosen test case
    
    switch prob_case
        case {'scalarMEG', 'vectorMEG'}
            numberMeasurements = device.MEG.numberMeasurement;
            scaleFac = 1e-16;
        case 'vectorEEG'
            numberMeasurements = device.EEG.numberMeasurement;
            scaleFac = 1;
    end
    
    switch data_generatingCase
        case 'const'
            % Constant test function F^* = = c € R
            caseParameter = {'c', 1e-6, 'caseType', 'const'};  %constant name c, constant value 1
        case 'onb'
            % ONB as test function F^* = G_(0,n,j)
            caseParameter = {'n',3,'j',2 , 'const', 1e-1,'caseType', 'onb'}; % spherical harmonics degree n and order j
        case '2splines'
            % Spline test function F^* = coeff_k* A^k K(*,y_k)
            spline_truncationIndex_data = 500; 
            h_data = 0.8; 
            n=1:spline_truncationIndex_data;
            sequence_data = n.*h_data.^(2*n);
            k=[10,min(80,numberMeasurements)];
            coeff = zeros(numberMeasurements,1);
            coeff(k) = -[7.5e-2,5e-2]*scaleFac; % yields realistic order of magnitude 
            caseParameter = {'k',k,'sequence', sequence_data, 'coeff', coeff, 'spline_truncationIndex_data', spline_truncationIndex_data, 'h', h_data, 'caseType', '2splines'};
        case '2splines_avoidIC'
            % Spline test function F^* = coeff_k* A^k K(*,y_k)
            spline_truncationIndex_data = 500; 
            h_data = 0.8; 
            n=1:spline_truncationIndex_data;
            sequence_data = h_data.^(2*n+2)./((2*n+5).*(2*n+1));
            k=[10,min(80,numberMeasurements)];
            coeff = zeros(numberMeasurements,1);
            coeff(k) =  -[7.5e-2,5e-2]*scaleFac;  % yields realistic order of magnitude 
            num_integration_point = 100;%4e6; % at least 100.000 points for acceptable numerical integration. more is better.
            caseParameter = {'k',k,'sequence', sequence_data, 'coeff', coeff, 'spline_truncationIndex_data', spline_truncationIndex_data,'num_integration_point', num_integration_point, 'h', h_data, 'integration points', num_integration_point, 'caseType', '2splines_avoidIC'};
        case '2splines_vec'
            spline_truncationIndex_data = 500; 
            h_data = 0.8; 
            n=1:spline_truncationIndex_data;
            sequence_data = n.*h_data.^(2*n);
            sequence_data = sequence_data.*(n.*(n+1))./(2*n+3).*(2*n+5);
            k=[10,min(80,numberMeasurements)];
            coeff = zeros(numberMeasurements,1);
            coeff(k) = -[6.5e3,4e3]*scaleFac; % yields realistic order of magnitude 
            caseParameter = {'k',k,'sequence', sequence_data, 'coeff', coeff, 'spline_truncationIndex_data', spline_truncationIndex_data, 'h', h_data, 'caseType', '2splines'};
  
    end
    
    % Calculate and plot the new exact solution of the test case
    generate_solution(device, data_generatingCase, caseParameter, prob_case, plot_program)
    if visualize_current && isequal(prob_case,'scalarMEG')
        generate_solution(device, data_generatingCase, caseParameter, 'vectorMEG', plot_program)
    end
    
    % Generate the data and noise it
    switch prob_case
        case {'scalarMEG', 'vectorMEG'}
            data = generate_data_scalarMEG(device, prob_case, data_generatingCase, caseParameter);
        case 'vectorEEG'
            data = generate_data_vectorEEG(device, data_generatingCase, caseParameter);
    end
end

function generate_solution(device, data_generatingCase, caseParameter, prob_case, plot_program)
    arguments
        device (1,1) struct
        data_generatingCase (1,:) char {mustBeMember(data_generatingCase,{'const','onb','2splines', '2splines_avoidIC', '2splines_vec'})} = 'const'
        caseParameter (1,:) cell  = {'c', 1, 'caseType', 'const'}
        prob_case (1,:) char {mustBeMember(prob_case, {'scalarMEG', 'vectorMEG','vectorEEG'})} = 'scalarMEG'
        plot_program (1,:) char {mustBeMember(plot_program, {'none','matlab', 'tikz'})} = 'tikz'
   end
    
    switch data_generatingCase
        case 'const'
            [eta_x,eta_y,eta_z,evalRadius,~,~] = generate_plotPositions(device);
            switch prob_case
                case 'scalarMEG'
                    % caseParameter = {'c', VAL};
                    valSolution = caseParameter{2}*ones(size(eta_x));
                case {'vectorMEG', 'vectorEEG'}
                    valSolution = caseParameter{2}*ones(3,size(eta_x,1),size(eta_x,2)); 
            end
        case 'onb'
            [eta_x,eta_y,eta_z,evalRadius,phi,t] = generate_plotPositions(device);
    
            % caseParameter = {'n', VAL, 'j', VAL}
            switch prob_case
                case 'scalarMEG'
                    spH = sphericalHarmonics(caseParameter{2}, phi, t);
                    valSolution = caseParameter{6}*sqrt((2*caseParameter{2}+5)/device.radiusCerebrum^3).*(evalRadius/device.radiusCerebrum).^(caseParameter{2}+1).*squeeze(spH(caseParameter{2}^2+caseParameter{2}+caseParameter{4}+1,:,:));
                case 'vectorMEG'
                    [eta_x,eta_y,eta_z,evalRadius,phi,t] = generate_plotPositions(device);
                    spherHarm = vecSpherHarm(3,caseParameter{2}, phi,t);
                    spherHarm = squeeze(spherHarm(caseParameter{2}^2+caseParameter{2}+caseParameter{4}+1,:,:,:));

                    valSolution = (2*sqrt((2*caseParameter{2}+5)*(2*caseParameter{2}+3)^2/(caseParameter{2}*(caseParameter{2}+1)))*caseParameter{6} * evalRadius^3 / device.radiusCerebrum^(11/2))*spherHarm;
                    valSolution = reshape(valSolution,3,size(eta_x,1),size(eta_x,2)); %fT conversion
                case 'vectorEEG'
                    valSolution = NaN*ones(3,size(eta_x,1),size(eta_x,2)); %case not implemented
            end
            case {'2splines', '2splines_avoidIC', '2splines_vec'}
            % In this case, plotPositions are already calculated in
            % evaluate_spline_XYZ function           
                switch prob_case
                    case 'scalarMEG'
                        % caseParameter = {'k', VAL,'sequence', VAL, 'coeff', VAL, 'spline_truncationIndex_data', VAL, 'h', VAL}
                        [eta_x,eta_y,eta_z,evalRadius,valSolution,~,~] = evaluate_spline(prob_case, device,caseParameter{4},caseParameter{6});
                    case {'vectorEEG','vectorMEG'}
                        switch prob_case
                            case 'vectorMEG'
                                % caseParameter = {'k', VAL,'sequence', VAL, 'coeff', VAL, 'spline_truncationIndex_data', VAL, 'h', VAL}
                                n = 1:length(caseParameter{4});
                                sequence = 4*(2*n+3).*(2*n+5).*caseParameter{4}./(device.radiusCerebrum^2.*(n+1).*n);
                            case 'vectorEEG'
                                sequence = caseParameter{4};
                        end
                        % caseParameter = {'k', VAL,'sequence', VAL, 'coeff', VAL, 'spline_truncationIndex_data', VAL, 'h', VAL}
                        [eta_x,eta_y,eta_z,evalRadius,valSolution_x, valSolution_y, valSolution_z] = evaluate_spline(prob_case, device,sequence,caseParameter{6});
                        valSolution = zeros(3,size(eta_x,1), size(eta_x,2));
                        valSolution(1,:,:) = reshape(valSolution_x, size(eta_x,1), size(eta_x,2));
                        valSolution(2,:,:) = reshape(valSolution_y, size(eta_x,1), size(eta_x,2));
                        valSolution(3,:,:) = reshape(valSolution_z, size(eta_x,1), size(eta_x,2));
                end
    end
    plotFile = struct('method', 'exactSolution', 'x', evalRadius*eta_x, 'y', evalRadius*eta_y, 'z', evalRadius*eta_z, 'valSolution', valSolution, 'data_generatingCase', data_generatingCase, 'prob_case', prob_case, 'caseParameter', {caseParameter});    
    save(['../solution/plot_' , data_generatingCase, '_', prob_case, '_exactSolution.mat'], '-struct', 'plotFile');
    
    switch prob_case
        case 'scalarMEG'
            plot_scalarFunction3d(eta_x*evalRadius,eta_y*evalRadius,eta_z*evalRadius, valSolution, ['../pics/', data_generatingCase,  '_', prob_case, '_exactSolution.tex'], plot_program)
        case {'vectorMEG', 'vectorEEG'}
            plot_scalarFunction2current3d(eta_x*evalRadius,eta_y*evalRadius,eta_z*evalRadius, valSolution, ['../pics/', data_generatingCase,  '_', prob_case, '_exactSolution.tex'], plot_program);
     end
end

function data = generate_data_vectorEEG(device, data_generatingCase, caseParameter)
    arguments
        device (1,1) struct
        data_generatingCase (1,:) char {mustBeMember(data_generatingCase,{'2splines'})} = '2splines'
        caseParameter (1,:) cell  = {'c', 1, 'caseType', 'const'}
    end

    splineMat = calculate_splineMatrix('vectorEEG', device, caseParameter{4});
    data = (caseParameter{6}(caseParameter{2})'*splineMat(caseParameter{2},:))'; %already in fT
     
    dataFile = struct('type', 'vectorEEG', 'U', data,  'UUnit', 'DETERMINE', 'No_s', length(data), 'Pos', device.EEG.measurementPosition', 'noiseLevel', 0, 'noiseType', 'none');
    save(['../solution/data_', data_generatingCase, '_vectorEEG_noiseLevel_0_noiseType_none.mat'], '-struct', 'dataFile'); 
end

function data = generate_data_scalarMEG(device, prob_case, data_generatingCase, caseParameter)
    arguments
        device (1,1) struct    
        prob_case (1,:) char {mustBeMember(prob_case, {'scalarMEG', 'vectorMEG', 'vectorEEG'})} = 'scalarMEG'
        data_generatingCase (1,:) char {mustBeMember(data_generatingCase,{'const','onb','2splines','2splines_avoidIC', '2splines_vec'})} = 'const'
        caseParameter (1,:) cell  = {'c', 1, 'caseType', 'const'}
    end
    % unit F^*_(n,j) = [F^*_(n,j)] = Am^(-1)
    % unit mu0 = [mu0] = NA^(-2)
    % unit data = [Bz] = NA^(-1)m^(-1) = T
            
    % Preallocation
    data = zeros(device.MEG.numberMeasurement,1);
    
    switch data_generatingCase
        case 'const'
            %Function has radial part F^*_(n,j)(r) = sqrt(4*pi) * const * delta_(n,0) * delta_(j,0)
            % y_(0,0) is in null space of forward operator, hence data equals zero
            disp('Data for constant test function calculated');
        case 'onb'
            % caseParameter = {'n', VAL, 'j', VAL}
            % Functional is an ONB with radial degree m=0 (null space), that is F = Q_0^(n+3/2)(R;r)Y_(n,j)(xi) = sqrt(2*n+5/R^3)(r/R)^(n+1) Y_(n,j)(xi)
            
            % radial function is F^*_(n,j)(r) = sqrt((2*n+5)/R^3)(r/R)^(n+1)
            % derivative of radial Function is F^*_(n,j)'(r)=sqrt((2*n+5)/R)(n+1)/R^2 (r/R)^n
            % A F^*= -mu0 * 1/sqrt((2n+1)(n+1))(F^*_(n,j)'(R)R - (n-1)F^*_(n,j)(R))(R/r)^(n+2) nu(r*xi)*y_(n,j)^(1)(xi)
            %      = -2 * mu0 sqrt((2n+5)/((2n+1)(n+1)R^3)) (R/r)^(n+2) nu(r*xi)*y_(n,j)^(1)(xi)
            
            coefficient = -2 * sqrt((2*caseParameter{2}+5)/((2*caseParameter{2}+1)*(caseParameter{2}+1))) * device.mu0 * device.radiusCerebrum^(-1.5);
            
            % Calculate Edmonds vector spherical harmonics of type 1 via linear combination of Morse-Feshbach vector spherical harmonics
            vecSpH1 = vecSpherHarm(1,caseParameter{2}, atan2(device.MEG.measurementPositionNormalized(:,2),device.MEG.measurementPositionNormalized(:,1)),device.MEG.measurementPositionNormalized(:,3));
            vecSpH1 = squeeze(vecSpH1(caseParameter{2}^2+caseParameter{4}+caseParameter{2}+1,:,:));
            vecSpH2 = vecSpherHarm(2,caseParameter{2}, atan2(device.MEG.measurementPositionNormalized(:,2),device.MEG.measurementPositionNormalized(:,1)),device.MEG.measurementPositionNormalized(:,3));
            vecSpH2 = squeeze(vecSpH2(caseParameter{2}^2+caseParameter{4}+caseParameter{2}+1,:,:));
            vecSpH = sqrt((caseParameter{2}+1)/(2*caseParameter{2}+1))*vecSpH1 - sqrt(caseParameter{2}/(2*caseParameter{2}+1))*vecSpH2;
            
            % Calculate the data
            data = caseParameter{6}*coefficient * (device.radiusCerebrum./device.MEG.measurementPositionRadius).^(caseParameter{2}+2).*dot(device.MEG.normals',vecSpH(:,:))';
            
            % Magnetic flux density of human brain in 50-500fT (femto T), transfere data to fT
            data = data*1e15;
        case {'2splines', '2splines_vec'}
            % Calculates the data values very fast by means of the svd, confesses inverse
            % crime (same method for creating data and solving inverse
            % problem used)
            splineMat = calculate_splineMatrix('scalarMEG', device, caseParameter{4},caseParameter{2});
            data = (caseParameter{6}(caseParameter{2})'*splineMat(caseParameter{2},:))'; %already in fT
            
         case '2splines_avoidIC'
            % Calculates the data via Monte Carlo Integration which avoids
            % the inverse crime but is significant slower and requires a
            % substaintial supplementary amount of memory.
            
            %splineMat = integragte_quasiMonteCarlo_ball_scalarMEG(device,caseParameter{10},caseParameter{14},caseParameter{14});
            %save(['../solution/splineMat', num2str(caseParameter{12}*100),'_AIC.mat'], 'splineMat');
            disp('load cluster result');
            load(['../solution/splineMat', num2str(caseParameter{12}*100),'_AIC.mat'], 'splineMat');
            data = (caseParameter{6}(caseParameter{2})'*splineMat(caseParameter{2},:))';
                    
    end
    dataFile = struct('type', prob_case, 'Bz', data,  'BzUnit', 'f NA^(-1)m^(-1)', 'No_s', length(data), 'Pos', device.MEG.measurementPosition', 'z_dir', device.MEG.normals', 'noiseLevel', 0, 'noiseType', 'none');
    save(['../solution/data_', data_generatingCase, '_', prob_case '_noiseLevel_0_noiseType_none.mat'], '-struct', 'dataFile'); 
end