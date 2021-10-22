function spline_MEEG(data_generatingCase, plot_program, h, spline_sequenceCase, spline_truncationIndex, verify_closedSeriesRepresentations, verify_basisMathFunctions, reg_method, rel_regPar, prob_case, recalc, noiseLevel,visualize_splineMatrix, visualize_current)
    % Code written by S. Leweke (Uni Siegen) started in 2020 based on
    % research code from dissertation in 2018. 
    % Details on MEEG modelling can be found
    % [1] Leweke, S. (2018). The Inverse Magneto-electroencephalography Problem for the Spherical Multiple-shell Model. Universität Siegen.
    % [2] Leweke, S., Michel, V., & Fokas, A. S. (2020). Electro-magnetoencephalography for a spherical multiple-shell model: Novel integral operators with singular-value decompositions. Inverse Probl., 36(3), 35003. https://doi.org/10.1088/1361-6420/ab291f
    % [3] Leweke, S., Michel, V. (2021). Vector-valued spline method for the spherical multiple-shell electro-magnetoencephalography problem.
    arguments
        data_generatingCase (1,:) char {mustBeMember(data_generatingCase, {'const','onb','2splines', '2splines_avoidIC', '2splines_vec','real'})} = '2splines'
        plot_program (1,:) char {mustBeMember(plot_program, {'none', 'matlab', 'tikz'})} = 'matlab'
        h (1,1) double {mustBeNumeric, mustBeReal, mustBeSmallerOne(h)} = 0.85
        spline_sequenceCase (1,:) char {mustBeMember(spline_sequenceCase,{'true_diss_spline_scalarMEG', 'stated_diss_spline_scalarMEG', 'sequence_scalarMEG', 'sequence_vectorMEG','sequence_vectorEEG'})} = 'sequence_vectorEEG'
        spline_truncationIndex (1,1) double {mustBeNumeric, mustBeReal, mustBePositive, mustBeInteger} = 200;
        verify_closedSeriesRepresentations (1,1) logical = false
        verify_basisMathFunctions (1,1) logical = false     
        reg_method (1,:) char {mustBeMember(reg_method, {'dont_solve','Tik', 'TV_PBP', 'none'})} = 'Tik'
        rel_regPar (1,:) double {mustBeReal} = logspace(-15,1,500)% 0-2% noise; logspace(-7,3,500) 5-10% noise 
        prob_case (1,:) char {mustBeMember(prob_case, {'scalarMEG', 'vectorMEG','vectorEEG'})} = 'scalarMEG'
        recalc (1,1) logical = false;
        noiseLevel (1,:) double {mustBeNumeric, mustBeReal, mustBeInteger} = 0;
        visualize_splineMatrix (1,1) logical  = false;
%                 The corresponding current can also be visualized with this feature. 
%                 Do only enable if necessary, since the function is comparable slow.
%                 The summation is the reason for the slowness.
        visualize_current (1,1) logical  = false; 
     end
    % Load device setting (sensor positions, radii and conductivity of tissues)
    close all;
    warning('off','all');
    device = load_device;
    
    % Correctness of subroutines can be verified, after changes in code:
    verify_subroutines(device, verify_closedSeriesRepresentations,verify_basisMathFunctions)
    
    if ~recalc && ~exist(['../solution/eval_splines_', prob_case, '.mat'],'file')
       recalc = true;
       disp('No splines evalable. Recalculation is required.')
    end
    
    if strcmp(prob_case,'vectorEEG') && ~(strcmp(data_generatingCase,'2splines') || strcmp(data_generatingCase,'real'))
        error('For vectorEEG only 2splines or real are available'),
    end
    % Build Spline Sequence
    % Series is truncataed at spline_truncationIndex and n=0 is not included in the series summation, due to null space structure of MEG and EEG problem.}
    % 'sequence' = kappan_n^(-2) from [1] Eq (20.11)
    % In [1] Sec 21.3.1. stated formulae for kappa_n containes a typo: case stated_diss_spline_scalarMEG is stated but true_diss_spline_scalarMEG is used in numerics

    disp('Initialize spline sequence');
    n=1:spline_truncationIndex;
    switch spline_sequenceCase
        case 'true_diss_spline_scalarMEG' 
            if h ~= 0.99
                warning('true spline_scalarMEG from [1] requires h=0.99')
            end
            sequence = h.^(n.*n).*n;
        case 'stated_diss_spline_scalarMEG' 
            if h ~= 0.99
                warning('stated_diss_spline_scalarMEG from [1] requires h=0.99')
            end
            sequence = h.^(n.*n)./n; 
        case 'sequence_scalarMEG'
            sequence = h.^n./n;
        case 'sequence_vectorEEG'
            sequence = h.^n;
        case 'sequence_vectorMEG'
            %sequence = h.^(n+3)./((n+500));
            sequence = h.^(6*n)./(n);
            sequence = sequence.*(2*n+3).*(2*n+5)./(n.*(n+1));
    end
   
   % Calculate the entire spline matrix
   splineMat = calculate_splineMatrix(prob_case, device, sequence);
   

    % Set of data at noisy pendands are calculated via chosen method
    switch data_generatingCase
        case 'real'
            disp('No data is generated. Load data instead.')           
        otherwise %'const','onb','2splines', '2splines_avoidIC', '2splines_vec'
            disp('Generate test case data and the exact solution')
            % Data is always generated based on the scalar problem. In the
            % vectorial case the corresponding exact solution current is calculated
            % additionally.
            data = generate_splineData(device, data_generatingCase, prob_case, plot_program, visualize_current);
            disp('Noise data')
            noise_data(data, prob_case, data_generatingCase,noiseLevel);
     end
              
    
    % Visualization of the spline matrix
    if visualize_splineMatrix
        switch plot_program
            case 'matlab'
                figure;
                imagesc(splineMat);
                colorbar;
                title('Visualization of Spline Matrix Entries.');
            case 'tikz'
                [X, Y] = meshgrid(linspace(1,102,102));
                export_splineMatrix2dat = [X(:), Y(:), splineMat(:)];
                save('../pics/splineMatrix.dat', 'export_splineMatrix2dat', '-ascii');
            case 'none'
                disp('Plotting is disabled.')
        end    
    end
    
%     Solve (regularized) inverse problem
    switch reg_method
        case {'Tik', 'TV_PBP', 'none'}
            for i=1:length(noiseLevel)
                % Load data depending on data case and noise level
                switch data_generatingCase
                    case {'const','onb','2splines', '2splines_avoidIC', '2splines_vec'} 
                        load_str = ['../solution/data_', data_generatingCase, '_', prob_case, '_noiseLevel_', num2str(noiseLevel(i)) ,'_noiseType_', 'gauss' ,'.mat'];
                    case 'real'
                        load_str = ['../data/data_', data_generatingCase, '_', prob_case,'.mat'];
                end
                if exist(load_str,'file')
                    dataFile = load(load_str);
                else
                    error('Expected data not available. Please check data input.')
                end
                % change if particular data set is supposed to be loaded
                % load_str =  'YORU PATH'
                % dataFile = load(load_str);
                switch prob_case
                    case {'scalarMEG', 'vectorMEG'}
                        data = dataFile.Bz;
                    case 'vectorEEG'
                        data = dataFile.U;
                end
                
                % Solve inverse spline problem
                disp('Solve inverse problem')
                list_regPar = max(max(abs(splineMat)))*rel_regPar;
                [coeff_reconst, list_regPar] = solve_inverseSplineProblem(data, splineMat, reg_method,list_regPar);
                approxFile = struct('type', reg_method, 'data', data, 'coeff_reconst', coeff_reconst, 'splineMat', splineMat, 'regPar',{list_regPar});
                save(['../solution/reconst_Tikhonov_', prob_case, '_data_generatingCase_', data_generatingCase, '_noise_', num2str(noiseLevel(i)), '_reg_Num_', num2str(length(list_regPar.regPar_all)), '_regParMin_', num2str(min(list_regPar.regPar_all)), '_regParMax_', num2str(max(list_regPar.regPar_all)),'.mat'], '-struct', 'approxFile'); 

                
                % Plot Spline
                disp('Plot best solutions')
                list_pcm = {'LCM auto', 'LCM man', 'DM', 'QOC', 'QCV'};
                if length(list_regPar) > 1
                    list_idx = [list_regPar.idx_auto,list_regPar.idx_man,list_regPar.idx_discrep,list_regPar.idx_QOC,list_regPar.idx_GCV];
                else
                    list_idx = 1;
                end
                coeff = coeff_reconst(:,list_idx);
                switch data_generatingCase
                    case 'real'
                        disp('Exact Solution does not exists.')
                        [~,idx,~] = plot_spline(device, sequence, coeff, 'approx', plot_program, data_generatingCase, prob_case, recalc, ['_h', num2str(h), '_','noise_', num2str(noiseLevel(i))]);
                        nrmse = NaN;
                    otherwise
                        [~,idx,nrmse] = plot_spline(device, sequence, coeff, 'diff', plot_program, data_generatingCase, prob_case, recalc, ['_h', num2str(h), '_','noise_', num2str(noiseLevel(i))]);
                        disp(['NRMSE: ', num2str(nrmse)]);
                end
                
                if visualize_current && isequal('scalarMEG', prob_case)
                    n = 1:length(sequence);
                    sequence = sequence.*(4*(2*n+3).*(2*n+5))./((n+1).*n.*device.radiusCerebrum^2);
                    [~,~,nrmse_vec] = plot_spline(device, sequence, coeff, 'diff', plot_program, data_generatingCase, 'vectorMEG', true, ['_current_h', num2str(h), '_','noise_', num2str(noiseLevel(i))]);
                    disp(['NRMSE vector: ', num2str(nrmse_vec)]);
                
                end
                
                % Save results
                filename = ['../solution/table_' , prob_case, '_noise', num2str(noiseLevel(i))];
                if exist([filename,'.mat'], 'file')
                    load(filename, 'A');
                else
                    load('../solution/table_empty.mat', 'A');
                end
                row_A = find(cell2mat(A(1,2:end)) == h, 1)+1;
                A{2,row_A} = norm((data-splineMat*coeff(:,idx)),2)/norm(data,2);
                A{5,row_A} = list_regPar.regPar_all(list_idx(idx))*1e-30;
                A{3,row_A} = nrmse;
                A{4,row_A} = list_pcm(idx);
                save(filename, 'A');
                disp(['Relative data misfit: ', num2str(A{2,row_A})]);
                disp(['Regularization parameter: ', num2str(A{5,row_A})]);
                disp(['idx: ', num2str(idx)])
            end
        case 'dont_solve'
            disp('Solving inverse problem is disabled.')
    end
end


function mustBeSmallerOne(a)
    if a >= 1
        eid = 'Val:notPermitted';
        msg = 'Spline parameter h must be smaller than one.';
        throwAsCaller(MException(eid,msg))
    end
end