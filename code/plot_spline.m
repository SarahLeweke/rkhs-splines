function [valApprox,idx, nrmse_idx] = plot_spline(device, sequence, coeff, method, plot_program, data_generatingCase, prob_case, recalc, nameTag)
    arguments
        device (1,1) struct
        sequence (1,:) double {mustBeNumeric,mustBeReal} 
        coeff (:,:) double {mustBeNumeric,mustBeReal} 
        method (1,:) char {mustBeMember(method,{'approx', 'diff'})} = 'approx'
        plot_program (1,:) char {mustBeMember(plot_program, {'none','matlab', 'tikz'})} = 'tikz'
        data_generatingCase (1,:) char {mustBeMember(data_generatingCase,{'none','const','onb','2splines','2splines_avoidIC', '2splines_vec', 'real'})} = 'const'
        prob_case (1,:) char {mustBeMember(prob_case, {'scalarMEG', 'vectorMEG','vectorEEG'})} = 'scalarMEG'
        recalc (1,1) logical = true;
        nameTag (1,:) char = '';
    end
    
    nrmse_idx = NaN;
    idx = 1;
    
    switch prob_case
         case 'scalarMEG'
            [eta_x,eta_y,eta_z, evalRadius, valApproxFull,~,~] = evaluate_spline(prob_case, device,sequence,coeff,recalc);
            if numel(valApproxFull) > numel(eta_x)
                minnorm = sqrt(sum((valApproxFull).^2,1)/size(valApproxFull,1));
                [nrmse_idx,idx] = min(minnorm);
                valApprox = reshape(valApproxFull(:,idx),size(eta_x,1),size(eta_x,2));
            else
                valApprox = valApproxFull;
                valApproxFull = valApproxFull(:);
            end
         case {'vectorMEG', 'vectorEEG'}
            [eta_x,eta_y,eta_z, evalRadius, valApproxFull_x, valApproxFull_y, valApproxFull_z] = evaluate_spline(prob_case, device,sequence,coeff,recalc);   
            valApprox_norm = sqrt(valApproxFull_x.^2 + valApproxFull_y.^2 + valApproxFull_z.^2);
            minnorm = sqrt(sum((valApprox_norm).^2,1)/size(valApprox_norm,1));
            [nrmse_idx,idx] = min(minnorm);
            valApprox = zeros(3,size(eta_x,1),size(eta_x,2));
            valApprox(1,:,:) = reshape(valApproxFull_x(:,idx),size(eta_x,1),size(eta_x,2));
            valApprox(2,:,:) = reshape(valApproxFull_y(:,idx),size(eta_x,1),size(eta_x,2));
            valApprox(3,:,:) = reshape(valApproxFull_z(:,idx),size(eta_x,1),size(eta_x,2));
    end
   
    switch method
        case 'diff'
            switch data_generatingCase
                case {'const','onb','2splines', '2splines_avoidIC', '2splines_vec'} 
                       load(['../solution/plot_', data_generatingCase, '_', prob_case, '_exactSolution.mat'],'valSolution')
            end
            
            switch prob_case
                case {'vectorMEG', 'vectorEEG'}
                    valSolution = reshape(valSolution,3,numel(eta_x));
                    norm_diff = sqrt((valApproxFull_x - valSolution(1,:)').^2 + (valApproxFull_y - valSolution(2,:)').^2 + (valApproxFull_z - valSolution(3,:)').^2);
                    
                    nrmse = sqrt(sum(norm_diff.^2,1)/numel(eta_x))/(max(max(sqrt(sum(valSolution.^2,1))))-min(min(sqrt(sum(valSolution.^2,1)))));
                    [nrmse_idx,idx] = min(nrmse);

                    valApprox = zeros(3,size(eta_x,1),size(eta_x,2));
                    valApprox(1,:,:) = reshape(valApproxFull_x(:,idx),size(eta_x,1),size(eta_x,2));
                    valApprox(2,:,:) = reshape(valApproxFull_y(:,idx),size(eta_x,1),size(eta_x,2));
                    valApprox(3,:,:) = reshape(valApproxFull_z(:,idx),size(eta_x,1),size(eta_x,2));

                    valSolution = reshape(valSolution,3,size(eta_x,1),size(eta_x,2));
                    valDiff = valApprox - valSolution;
                    
                    fileName = ['../pics/plot_' data_generatingCase '_' prob_case, nameTag, '_diff.tex'];
                    plot_scalarFunction2current3d(eta_x*evalRadius,eta_y*evalRadius,eta_z*evalRadius, valDiff, fileName, plot_program);
                    
                case 'scalarMEG'
                    nrmse = sqrt(sum((valApproxFull-valSolution(:)).^2,1)/numel(valSolution))/(max(max(abs(valSolution)))-min(min(abs(valSolution))));
                    [nrmse_idx,idx] = min(nrmse);

                    valApprox = reshape(valApproxFull(:,idx),size(valSolution,1),size(valSolution,2));
                    valDiff = abs(valApprox-valSolution);

                    fileName = ['../pics/plot_' data_generatingCase '_' prob_case, nameTag, '_diff.tex'];
                    plot_scalarFunction3d(eta_x*evalRadius,eta_y*evalRadius,eta_z*evalRadius, valDiff, fileName, plot_program);    
            end           
   end
   fileName = ['../pics/plot_' data_generatingCase '_' prob_case,  nameTag,'_approx.tex'];
   switch prob_case
       case {'vectorMEG','vectorEEG'}
           plot_scalarFunction2current3d(eta_x*evalRadius,eta_y*evalRadius,eta_z*evalRadius, valApprox, fileName, plot_program);
       case 'scalarMEG'
           plot_scalarFunction3d(eta_x*evalRadius,eta_y*evalRadius,eta_z*evalRadius, valApprox, fileName, plot_program);
   end
       
end