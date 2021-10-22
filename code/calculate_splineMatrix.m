function splineMat = calculate_splineMatrix(prob_case, device, sequence, list_column)
    % case list_column = NaN: calculates only left under triangle matrix
    % fist, than mirrors right-upper matrix. Fastes method if full matrix
    % required
    
    % case list_column = [10, 36] (e.g.): calculates only 10th and 36th
    % column of matrix. Faster if only up to 75% of the columns are
    % required
    arguments
        prob_case (1,:) char {mustBeMember(prob_case, {'scalarMEG', 'vectorMEG','vectorEEG'})} 
        device (1,1) struct
        sequence (1,:) double {mustBeNumeric, mustBeReal} 
        list_column (:,1) double {mustBeNumeric, mustBeReal} = NaN
    end
    
    switch prob_case
        case 'scalarMEG'
            % spline matrix is calculated by applying the forward functionals onto
            % the reproducing kernel twice, see Eq. (20.11) from [1] with
            % abbreviation 'sequence' = kappan_n^(-2) 

            % To speed up summation, the Edmonds vector spherical harmonic terms
            % are replaced by Morse-Feshbach vector spherical harmonics, addition
            % theorems and 3d-orthogonality considerations are taken into account.
            % This result is stated in [3] sec. 3

            % Preallocation
            n=1:length(sequence);
            splineMat = zeros(device.MEG.numberMeasurement);

            % nu(y_k)*(eta_k)
            norcomp = dot(device.MEG.normals,device.MEG.measurementPositionNormalized,2);

            if isnan(list_column)
                % For sake of speedup the inner for-loop is replaced by vectorization
                for j=1:device.MEG.numberMeasurement
                    % Spline matrix is symmetric, hence only triangle matrix is calculated 
                    k = 1:j;
                    splineMat(j,k) = calculate_splineColumn_scalarMEG(device, n, sequence, norcomp,j,k);
                end  

                % Calculate spline matrix values in fT instead of T ( functional is applied twice)
                splineMat = (splineMat + splineMat' - diag(diag(splineMat)))*device.mu0^2/(pi*device.radiusCerebrum^3)*1e30;
                save('../solution/splineMat.mat', 'splineMat');
            else
                for i=1:length(list_column)
                    j = list_column(i);
                    % Spline matrix is symmetric, hence only triangle matrix is calculated 
                    k = 1:device.MEG.numberMeasurement;

                    splineMat(j,k) = calculate_splineColumn_scalarMEG(device, n, sequence, norcomp,j,k);
                end  
                % Calculate spline matrix values in fT instead of T ( functional is applied twice)
                splineMat = splineMat*device.mu0^2/(pi*device.radiusCerebrum^3)*1e30;
                % only full spline matrix is saved to avoid confusion
            end
        case 'vectorMEG'
            % spline matrix is calculated by applying the forward functionals onto
            % the reproducing kernel twice, see Eq. (20.11) from [1] with
            % abbreviation 'sequence' = kappan_n^(-2) 

            % To speed up summation, the Edmonds vector spherical harmonic terms
            % are replaced by Morse-Feshbach vector spherical harmonics, addition
            % theorems and 3d-orthogonality considerations are taken into account.
            % This result is stated in [3] sec. 3

            % Preallocation
            n=1:length(sequence);
            splineMat = zeros(device.MEG.numberMeasurement);

            % nu(y_k)*(eta_k)
            norcomp = dot(device.MEG.normals,device.MEG.measurementPositionNormalized,2);

            for j=1:device.MEG.numberMeasurement
                % Spline matrix is symmetric, hence only triangle matrix is calculated 
                k = 1:device.MEG.numberMeasurement;

                splineMat(j,k) = calculate_splineColumn_vectorMEG(device, n, sequence, norcomp,j,k);
            end  
            % Calculate spline matrix values in fT instead of T ( functional is applied twice)
            splineMat = splineMat*device.mu0^2/(4*pi*device.radiusCerebrum)*1e30;
        case 'vectorEEG'
            % spline matrix is calculated by applying the forward functionals onto
            % the reproducing kernel twice, see Sec. 20.3 from [1] with
            % abbreviation 'sequence' = kappan_n^(-2) 

            % To speed up summation, the Edmonds vector spherical harmonic terms
            % are replaced by Morse-Feshbach vector spherical harmonics, addition
            % theorems and 3d-orthogonality considerations are taken into account.
            % This result is stated in [3] Proof Thm. 4.3

            % Preallocation
            splineMat = zeros(device.EEG.numberMeasurement);
            beta = calculate_beta(device, length(sequence), 'resive');
            sequence = sequence(1:length(beta));
            coefClen = zeros(length(beta)+1,1);
            n = 1:length(beta);

            % Load the measurement positions and calculate their radii
            measurementPosition = device.EEG.measurementPosition';
            measureRadius = sqrt(sum(measurementPosition.^2,1));
            measureNormalized = (measurementPosition./measureRadius)';

            for j=1:length(measurementPosition)
                for k=1:j
                    coefClen(2:end) = sequence.*(device.radiusCerebrum^2./(measureRadius(k)*measureRadius(j))).^(n+1).*((n+1).*(measureRadius(k)/device.radiusScalp).^(2*n+1)+ n).*...
                        ((n+1).*(measureRadius(j)/device.radiusScalp).^(2*n+1)+ n).*beta.*beta.*(2*n+1)./n;
                    splineMat(j,k) = clenshawLegendre(coefClen, measureNormalized(j,:)*measureNormalized(k,:)')/(4*pi*device.radiusCerebrum);
                    splineMat(k,j) = splineMat(j,k);
                end
            end
    end
            
end

function column = calculate_splineColumn_scalarMEG(device, n, sequence, norcomp,j,k)
    % Coefficients and points for Clenshaw algorithm evaluation are
    % calculated
    coefClen = zeros(length(k),length(sequence)+1);
    coefClen(:,2:end) = sequence.*(device.radiusCerebrum^2./(device.MEG.measurementPositionRadius(k)*device.MEG.measurementPositionRadius(j))).^(n+2).*(2*n+5)./((n+1).*(2*n+1));
    t = device.MEG.measurementPositionNormalized(j,:)*device.MEG.measurementPositionNormalized(k,:)';

    % Calculate all Clenshaw summations
    term1 = clenshawLegendre((coefClen.*[0,n+1])', t')'.*norcomp(k)*norcomp(j);
    sum_LegendreDiff = clenshawLegendreDiff(coefClen', t')';
    term2 = -norcomp(k).*sum_LegendreDiff.*(device.MEG.normals(j,:)*device.MEG.measurementPositionNormalized(k,:)'- t*norcomp(j))';
    term3 = -norcomp(j)*sum_LegendreDiff.*(device.MEG.normals(k,:)*device.MEG.measurementPositionNormalized(j,:)'- t'.*norcomp(k));
    sum_LegendreDiff2 = clenshawLegendreDiff((coefClen.*[0,1./(n+1)])', t')'; %Die sind korrekt
    sum_LegendreSecondDiff = clenshawLegendreDiff2((coefClen.*[0,0,1./(n(2:end)+1)])', t')'; %Die sind korrekt

    % Collect all summands
    term5 = (sum_LegendreDiff2 + t'.*sum_LegendreSecondDiff).*(norcomp(k).*norcomp(j).*t' ...
            - norcomp(k).*(device.MEG.normals(j,:)*device.MEG.measurementPositionNormalized(k,:)')' ...
            - (device.MEG.normals(k,:)*device.MEG.measurementPositionNormalized(j,:)').*norcomp(j))...
        + (device.MEG.normals(k,:)*device.MEG.measurementPositionNormalized(j,:)').*(device.MEG.normals(j,:)*device.MEG.measurementPositionNormalized(k,:)')'.*sum_LegendreSecondDiff ...
        + (device.MEG.normals(j,:)*device.MEG.normals(k,:)')'.*sum_LegendreDiff2;
  
    column = (term1+term2+term3+term5);
end

function column = calculate_splineColumn_vectorMEG(device, n, sequence, norcomp,j,k)
    % Coefficients and points for Clenshaw algorithm evaluation are
    % calculated
    coefClen = zeros(length(k),length(sequence)+1);
    coefClen(:,2:end) = sequence.*(device.radiusCerebrum^2./(device.MEG.measurementPositionRadius(k)*device.MEG.measurementPositionRadius(j))).^(n+2).*(n)./((2*n+3).*(2*n+1));
    t = device.MEG.measurementPositionNormalized(j,:)*device.MEG.measurementPositionNormalized(k,:)';

    % Calculate all Clenshaw summations
    term1 = clenshawLegendre((coefClen.*[0,n+1])', t')'.*norcomp(k)*norcomp(j);
    sum_LegendreDiff = clenshawLegendreDiff(coefClen', t')';
    term2 = -norcomp(k).*sum_LegendreDiff.*(device.MEG.normals(j,:)*device.MEG.measurementPositionNormalized(k,:)'- t*norcomp(j))';
    term3 = -norcomp(j)*sum_LegendreDiff.*(device.MEG.normals(k,:)*device.MEG.measurementPositionNormalized(j,:)'- t'.*norcomp(k));
    sum_LegendreDiff2 = clenshawLegendreDiff((coefClen.*[0,1./(n+1)])', t')'; 
    sum_LegendreSecondDiff = clenshawLegendreDiff2((coefClen.*[0,0,1./(n(2:end)+1)])', t')'; 

    % Collect all summands
    term5 = (sum_LegendreDiff2 + t'.*sum_LegendreSecondDiff).*(norcomp(k).*norcomp(j).*t' ...
            - norcomp(k).*(device.MEG.normals(j,:)*device.MEG.measurementPositionNormalized(k,:)')' ...
            - (device.MEG.normals(k,:)*device.MEG.measurementPositionNormalized(j,:)').*norcomp(j))...
        + (device.MEG.normals(k,:)*device.MEG.measurementPositionNormalized(j,:)').*(device.MEG.normals(j,:)*device.MEG.measurementPositionNormalized(k,:)')'.*sum_LegendreSecondDiff ...
        + (device.MEG.normals(j,:)*device.MEG.normals(k,:)')'.*sum_LegendreDiff2;
  
    column = (term1+term2+term3+term5);
end