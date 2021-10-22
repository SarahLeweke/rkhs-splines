function [eta_x,eta_y,eta_z, evalRadius, valApprox_x, valApprox_y, valApprox_z] = evaluate_spline(prob_case, device,sequence,coeff, recalc)
    arguments
        prob_case (1,:) char {mustBeMember(prob_case, {'scalarMEG', 'vectorMEG','vectorEEG'})} 
        device (1,1) struct
        sequence (1,:) double {mustBeNumeric,mustBeReal} 
        coeff (:,:) double {mustBeNumeric,mustBeReal, mustBeDimFit(coeff,device,prob_case)}
        recalc (1,1) logical = true;
    end
    
    %Load positions for plotting
    [eta_x,eta_y,eta_z,evalRadius,~,~] = generate_plotPositions(device);

    % Loads the calculated splines, if they are available 
    if recalc
        splines = calculate_splines(prob_case, device, sequence, coeff, eta_x,eta_y,eta_z,evalRadius);

        evalSplineFile = struct('splines', splines, 'x', evalRadius*eta_x, 'sequence', sequence);
        save(['../solution/eval_splines_', prob_case, '.mat'], '-struct', 'evalSplineFile'); 
    else
        evalSplines = load(['../solution/eval_splines_', prob_case, '.mat']); 
        if isequal(evalSplines.x, eta_x*evalRadius) && isequal(evalSplines.sequence,sequence)
            splines = evalSplines.splines;
        else
            disp('splines do not match. need to recalc.')
            splines = calculate_splines(prob_case, device, sequence, coeff, eta_x,eta_y,eta_z,evalRadius);
            evalSplineFile = struct('splines', splines, 'x', evalRadius*eta_x, 'sequence', sequence);
            save(['../solution/eval_splines_', prob_case, '.mat'], '-struct', 'evalSplineFile'); 
        end      
    end

    switch prob_case
        case 'scalarMEG'
            % Multiplication with the spline coefficients.
            if size(coeff,2) == 1
                valApprox = reshape(coeff'*splines, size(eta_x,1), size(eta_x,2));
            else
                valApprox = zeros(numel(eta_x),size(coeff,2));
                for j=1:size(coeff,2)
                    valApprox(:,j) = coeff(:,j)'*splines;
                end
            end
            valApprox_x = valApprox;
            valApprox_y = [];
            valApprox_z = [];
        case {'vectorMEG', 'vectorEEG'}
            % Multiplication with the spline coefficients.
            valApprox_x = zeros(numel(eta_x),size(coeff,2));
            valApprox_y = zeros(numel(eta_x),size(coeff,2));
            valApprox_z = zeros(numel(eta_x),size(coeff,2));
            for j=1:size(coeff,2)
                valApprox_x(:,j) = coeff(:,j)'*squeeze(splines(:,1,:));
                valApprox_y(:,j) = coeff(:,j)'*squeeze(splines(:,2,:));
                valApprox_z(:,j) = coeff(:,j)'*squeeze(splines(:,3,:));
            end
    end
end

function splines = calculate_splines(prob_case, device, sequence, coeff, eta_x,eta_y,eta_z,evalRadius) %scalarMEG
    % Preallocation
    switch prob_case
        case 'scalarMEG'
            n=1:length(sequence);
            coefClen = zeros(1,length(sequence)+1);
            val1 = zeros(device.MEG.numberMeasurement,length(eta_x(:))); 
            val2 = zeros(device.MEG.numberMeasurement,length(eta_x(:))); 

            if size(coeff,2) == 1
                list_index = find(coeff ~=0);
            else
                list_index = 1:device.MEG.numberMeasurement;
            end

            % Calculate only the splines who influence the solution. speedup in the
            % '2splines' test case is approx 37.
            for k=1:length(list_index)
                t = [eta_x(:),eta_y(:),eta_z(:)]*device.MEG.measurementPositionNormalized(list_index(k),:)';
                coefClen(2:end) = sequence.*(2*n+5).*(evalRadius/device.MEG.measurementPositionRadius(list_index(k))).^(n+1)/device.MEG.measurementPositionRadius(list_index(k));
                val1(list_index(k),:) = clenshawLegendre(coefClen', t').*(device.MEG.normals(list_index(k),:)*device.MEG.measurementPositionNormalized(list_index(k),:)');
                val2(list_index(k),:) = - clenshawLegendreDiff(coefClen'.*[0,1./(n+1)]', t').*(device.MEG.normals(list_index(k),:)*[eta_x(:)';eta_y(:)';eta_z(:)'] - t'.*(device.MEG.normals(list_index(k),:)*device.MEG.measurementPositionNormalized(list_index(k),:)'));
            end
            % Calculate Spline in fT instead of T
            splines = -(device.mu0*(val1+val2)./(2*pi*device.radiusCerebrum^2))*10^15;   
        case 'vectorMEG'
            % Preallocation
            n=1:length(sequence);
            coefClen = zeros(1,length(sequence)+1);
            val1 = zeros(device.MEG.numberMeasurement,3,length(eta_x(:))); 
            val2 = zeros(device.MEG.numberMeasurement,3,length(eta_x(:))); 
            val3 = zeros(device.MEG.numberMeasurement,3,length(eta_x(:))); 

            if size(coeff,2) == 1
                list_index = find(coeff ~=0);
            else
                list_index = 1:device.MEG.numberMeasurement;
            end

            % Calculate only the splines who influence the solution. speedup in the
            % '2splines' test case is approx 37.
            for k=1:length(list_index)
                t = [eta_x(:),eta_y(:),eta_z(:)]*device.MEG.measurementPositionNormalized(list_index(k),:)';
                coefClen(2:end) = sequence.*(evalRadius/device.MEG.measurementPositionRadius(list_index(k))).^n/(device.MEG.measurementPositionRadius(list_index(k))^2);
                nuDotEta = device.MEG.normals(list_index(k),:)*device.MEG.measurementPositionNormalized(list_index(k),:)';
                nuProject = device.MEG.normals(list_index(k),:)' - nuDotEta.*device.MEG.measurementPositionNormalized(list_index(k),:)';
                val1(list_index(k),:,:) = clenshawLegendreDiff(coefClen', t').* nuDotEta .* cross([eta_x(:)';eta_y(:)';eta_z(:)'], repmat(device.MEG.measurementPositionNormalized(list_index(k),:)',1,numel(eta_x)),1);
                val2(list_index(k),:,:) = - clenshawLegendreDiff(coefClen'.*[0,1./(n+1)]', t').*cross([eta_x(:)';eta_y(:)';eta_z(:)'], repmat(nuProject,1,numel(eta_x)),1);
                val3(list_index(k),:,:) = - clenshawLegendreDiff2(coefClen'.*[0,1./(n+1)]', t').* ([eta_x(:)';eta_y(:)';eta_z(:)']'*nuProject)'.*cross([eta_x(:)';eta_y(:)';eta_z(:)'], repmat(device.MEG.measurementPositionNormalized(list_index(k),:)',1,numel(eta_x)),1);%.*(nuProject.*device.MEG.measurementPositionNormalized(list_index(k),:)');
            end
            % Calculate Spline in fT instead of T
            splines = -(device.mu0*(val1+val2 + val3)./(4*pi))*10^15;  
        case 'vectorEEG'
             if size(coeff,2) == 1
                list_index = find(coeff ~=0);
            else
                list_index = 1:device.EEG.numberMeasurement;
            end

            % Load the measurement positions and calculate their radii
            measurementPosition = device.EEG.measurementPosition;
            measureRadius = sqrt(sum(measurementPosition.^2,2));
            measureNormalized = (measurementPosition./measureRadius);

            % Preallocation
            beta = calculate_beta(device, length(sequence), 'resive');
            sequence = sequence(1:length(beta));
            n=1:length(sequence);
            coefClen = zeros(1,length(sequence)+1);
            splines = zeros(length(measureRadius),3,numel(eta_x)); 


            for k=1:length(list_index)
                t = [eta_x(:),eta_y(:),eta_z(:)]*measureNormalized(list_index(k),:)';
                coefClen(2:end) = sequence.*(2*n+1).*beta.*((n+1).*(measureRadius(list_index(k))/device.radiusScalp).^(2*n+1) + n).*(evalRadius/measureRadius(list_index(k))).^(n-1)/measureRadius(list_index(k))^2;
                splines(list_index(k),:,:) = (clenshawLegendre(coefClen', t)'.*[eta_x(:),eta_y(:),eta_z(:)] + clenshawLegendreDiff(coefClen'.*[0;1./n'], t)'.*(measureNormalized(list_index(k),:) - t.*[eta_x(:),eta_y(:),eta_z(:)]))';
            end
            splines = reshape(splines./(4*pi), length(measureRadius), 3, size(eta_x,1), size(eta_y,2));
    end
end
    

function mustBeDimFit(a,b,c)
    switch c 
        case {'scalarMEG','vectorMEG'}
            if size(a,1) ~= b.MEG.numberMeasurement
                eid = 'Size:notPermitted';
                msg = 'Data and coefficient dimension must fit.';
                throwAsCaller(MException(eid,msg))
            end
        case 'vectorEEG'
            if size(a,1) ~= b.EEG.numberMeasurement
                eid = 'Size:notPermitted';
                msg = 'Data and coefficient dimension must fit.';
                throwAsCaller(MException(eid,msg))
            end            
    end
end