function [regPar_best,x_best,idx] = finding_best_regularizationParameter(method, regPar, data, A, x, par)
    arguments
        method (1,:) char {mustBeMember(method,{'lCurve_auto','lCurve_man','discrepancy','QOC','GCV'})}
        regPar (1,:) double  {mustBeNumeric,mustBeReal, mustBePositive} 
        data (:,1) double  {mustBeNumeric,mustBeReal} 
        A (:,:) double  {mustBeNumeric,mustBeReal} 
        x (:,:) double  {mustBeNumeric,mustBeReal, mustBeSizeFit(data, A, x, regPar)} 
        par (1,4) cell = {'N2S', 0.000001, 'tau', 1.05}
    end
    
    % Sorting regParameter descending
    [regPar, I] = sort(regPar, 'descend');
    for i=1:size(x,1)
        x(i,:) = x(i,I);
    end
    
    % Calculate data misfits 
    data_error_rel = sqrt(sum(((A*x-data)./data).^2,1));
    data_error = sqrt(sum((A*x-data).^2,1));
    recon_norm = sqrt(sum(x.^2,1));
                    
    switch method
        case 'lCurve_auto'
            [~,idx] = min(sqrt(sum((A*x-data).^2,1)).*sqrt(sum(x.^2,1)));
        case 'lCurve_man'
            %loglog(recon_norm,data_error,'.');
            figure;
            loglog(recon_norm,data_error_rel,'o');
            labels = 1:length(regPar);
            labelpoints(recon_norm,data_error_rel,labels,'SW',0.2,1)
            prompt = {'Enter the index of the LCurve knick '};
            dlgtitle = 'Index';
            definput = {num2str(floor(length(regPar)/2))};
            dims = [1 length(regPar)];
            disp('Please choose the knick point. If you are ready for input, press enter')
            pause;
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            idx = str2double(answer{1});
        case 'discrepancy'
            noiseLvl = par{2}*norm(data)/sqrt(length(data));
            idx = find(data_error <= par{4}*noiseLvl*sqrt(length(data)),1,'first');
        case 'QOC'
            singVal = svd(A);
            rho = sqrt(sum((singVal*ones(size(regPar))./(singVal.^2*ones(size(regPar)) + ones(size(singVal))*regPar)).^2,1));
            [~,idx_max] = max(rho);
            K_max = find(rho < 0.5*rho(idx_max),1,'last');
            if isempty(K_max)
                K_max = length(rho)-1;
            end
            [~,idx] = min(sum((x(:,1:K_max)-x(:,2:(K_max+1))).^2,1));
        case 'GCV'            
            singVal = svd(A);
            rho = sqrt(sum((singVal*ones(size(regPar))./(singVal.^2*ones(size(regPar)) + ones(size(singVal))*regPar)).^2,1));
            [~,idx_max] = max(rho);
            K_max = find(rho < 0.5*rho(idx_max),1,'last');
            if isempty(K_max)
                K_max = length(rho)-1;
            end
            tmp = zeros(1,K_max);
            for k=1:K_max
               R = (A'*A+diag(regPar(k)))^(-1)*A';
               tmp(k) = norm(A*x(:,k)-data)/(trace(eye(length(data))-A*R)./length(data)).^2;
            end
            [~,idx] = min(tmp);
    end
    x_best = x(:,idx);
    regPar_best = regPar(idx);
    
    % Retourns index with the original sorting
    lst = 1:length(regPar);
    idx = lst(I(idx));
end

function mustBeSizeFit(g,A,x, regPar)
% Test for fitting dimension
    if ~isequal(size(A,1), size(A,2)) || ~isequal(size(A,2), size(x,1)) || ~isequal(size(A,1), length(g))
        eid = 'Size:notEqual';
        msg = 'Dimensions of Spline Problem related variables does not match.';
        throwAsCaller(MException(eid,msg))
    end
    if ~isequal(size(x,2),length(regPar))
        eid = 'Size:notEqual';
        msg = 'Number of regularization parameters does not match size of reconstruction.';
        throwAsCaller(MException(eid,msg))
    end
end