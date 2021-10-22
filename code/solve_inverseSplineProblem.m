function [coeff_reconst, regPar_out] = solve_inverseSplineProblem(data, splineMat, reg_method, list_regPar)
    arguments
       data (:,1) double {mustBeReal}
       splineMat (:,:) double {mustBeReal, mustBeMatchData(splineMat, data)} 
       reg_method (1,:) char {mustBeMember(reg_method, {'Tik', 'TV_PBP', 'none'})} = 'Tik'
       list_regPar (1,:) double {mustBeReal} = max(max(abs(splineMat)))*linspace(3e-2,0.5,500); %logspace(-15,1,500)% 0-2% noise; logspace(-7,3,500) 5-10% noise
    end
   
    regPar_out = struct('regPar_all', list_regPar);
    switch reg_method
        case 'none'
            coeff_reconst = linsolve(splineMat, data);
        case 'Tik'
            % Spline Problem = arg min (||g-SplineMat*x||^2 + lambda ||x||^1)
            coeff_reconst = zeros(size(splineMat,1),length(list_regPar));
            for i = 1:length(list_regPar)
                coeff_reconst(:,i) = solveTikhonovReg(data,splineMat,list_regPar(i));
            end
            labels = {};
            if length(list_regPar) > 5
                [~,~,regPar_out.idx_auto] = finding_best_regularizationParameter('lCurve_auto', list_regPar, data, splineMat, coeff_reconst);
                if ~isempty(regPar_out.idx_auto)
                    labels{length(labels)+1} = 'auto';
                end
                %[~,~,regPar_out.idx_man] = finding_best_regularizationParameter('lCurve_man', list_regPar, data, splineMat, coeff_reconst);
                regPar_out.idx_man = []; 
                if ~isempty(regPar_out.idx_man)
                    labels{length(labels)+1} = 'man';
                end
                [~,~,regPar_out.idx_discrep] = finding_best_regularizationParameter('discrepancy', list_regPar, data, splineMat, coeff_reconst);
                if ~isempty(regPar_out.idx_discrep)
                    labels{length(labels)+1} = 'discrep';
                end
                [~,~,regPar_out.idx_QOC] = finding_best_regularizationParameter('QOC', list_regPar, data, splineMat, coeff_reconst);
                if ~isempty(regPar_out.idx_QOC)
                    labels{length(labels)+1} = 'QOC';
                end
                [~,~,regPar_out.idx_GCV] = finding_best_regularizationParameter('GCV', list_regPar, data, splineMat, coeff_reconst);
                if ~isempty(regPar_out.idx_GCV)
                    labels{length(labels)+1} = 'idx_GCV';
                end
                data_error_rel = sqrt(sum(((splineMat*coeff_reconst-data)./data).^2,1));
                recon_norm = sqrt(sum(coeff_reconst.^2,1));
                figure;
                loglog(recon_norm,data_error_rel,'o');
                labelpoints(recon_norm([regPar_out.idx_auto,regPar_out.idx_man,regPar_out.idx_discrep,regPar_out.idx_QOC,regPar_out.idx_GCV]),data_error_rel([regPar_out.idx_auto,regPar_out.idx_man,regPar_out.idx_discrep,regPar_out.idx_QOC,regPar_out.idx_GCV]),labels,'SW',0.2,1)
            end
        case 'TV_PBP'
             % Solves  = arg min (0.5*||g-SplineMat*x||^2 + lambda ||x||)
            coeff_reconst = zeros(size(splineMat,1),length(list_regPar));
            for i = 1:length(list_regPar)
                coeff_reconst(:,i) = solveTVReg_penalizedBasisPursuit(data',splineMat,list_regPar(i));  
            end
    end
end

function x = solveTikhonovReg(b,A,lambda)
    x = (A'*A + lambda^2*eye(size(A'*A)))\(A'*b);
end

function mustBeMatchData(a,b)
    if size(a,1) ~= size(a,2) || size(a,1) ~= length(b)
        eid = 'Size:notPermitted';
        msg = 'Size of spline matrix and data must match.';
        throwAsCaller(MException(eid,msg))
    end
end