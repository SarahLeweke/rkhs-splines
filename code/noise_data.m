function noise_data(data, prob_case, data_generatingCase, noiseLevel ,noiseType)
    arguments
        data (:,1) double {mustBeNumeric, mustBeReal}   
        prob_case (1,:) char {mustBeMember(prob_case, {'scalarMEG', 'vectorMEG', 'vectorEEG'})} = 'scalarMEG'
        data_generatingCase (1,:) char {mustBeMember(data_generatingCase,{'none','const','onb','2splines','2splines_avoidIC', '2splines_vec'})} = 'const'
        noiseLevel (1,:) double {mustBeNumeric,mustBeReal,mustBePercent(noiseLevel)} = 1    
        noiseType (1,:) char {mustBeMember(noiseType,{'gauss'})} = 'gauss'
    end
    dataNoisy = generate_noisyData(data, noiseLevel, noiseType);
    save_noisyData(dataNoisy, prob_case,noiseLevel,noiseType,data_generatingCase);
end

function dataNoisy = generate_noisyData(data, noiseLevel, noiseType)
    arguments
        data (:,1) double {mustBeNumeric,mustBeReal}
        noiseLevel (1,:) double {mustBeNumeric,mustBeReal} = 1
        noiseType (1,:) char {mustBeMember(noiseType,{'gauss'})} = 'gauss'
    end
    
    dataNoisy = zeros(length(data), length(noiseLevel));
    switch noiseType
        case 'gauss'
            for i= 1:length(noiseLevel)
                if noiseLevel(i) > 0
                    % additive Gaussian white noise, mean 0, std depending
                    % relatively on noise level
                    dataNoisy(:,i) = data + normrnd(0,abs(data)*noiseLevel(i)/100,length(data),1);
                else
                    dataNoisy(:,i) = data;
                end
            end
    end
end

function save_noisyData(dataNoisy,prob_case,noiseLevel,noiseType,data_generatingCase)
    arguments
        dataNoisy double {mustBeNumeric,mustBeReal}
        prob_case (1,:) char {mustBeMember(prob_case, {'scalarMEG', 'vectorMEG', 'vectorEEG'})} = 'scalarMEG'
        noiseLevel (1,:) double {mustBeNumeric,mustBeReal,mustBeDimensionFit(dataNoisy,noiseLevel)} = 1
        noiseType (1,:) char {mustBeMember(noiseType,{'gauss'})} = 'gauss'
        data_generatingCase (1,:) char {mustBeMember(data_generatingCase,{'none','const','onb','2splines','2splines_avoidIC','2splines_vec'})}  = 'const'
    end
    
    dataFile = load(['../solution/data_', data_generatingCase, '_' , prob_case, '_noiseLevel_0_noiseType_none.mat']);
    
    switch prob_case
        case {'scalarMEG', 'vectorMEG'}
            for i=1:length(noiseLevel)
                dataFile.Bz = dataNoisy(:,i);
                dataFile.noiseLevel = noiseLevel(i);
                dataFile.noiseType = noiseType;
                save(['../solution/data_', data_generatingCase, '_' , prob_case, '_noiseLevel_' , num2str(noiseLevel(i)), '_noiseType_', noiseType ,'.mat'], '-struct', 'dataFile'); 
            end
        case 'vectorEEG'
            for i=1:length(noiseLevel)
                dataFile.U = dataNoisy(:,i);
                dataFile.noiseLevel = noiseLevel(i);
                dataFile.noiseType = noiseType;
                save(['../solution/data_', data_generatingCase, '_' , prob_case, '_noiseLevel_' , num2str(noiseLevel(i)), '_noiseType_', noiseType ,'.mat'], '-struct', 'dataFile'); 
            end    
    end
end

function mustBeDimensionFit(a,b)
% Test for fitting dimension
    if ~isequal(size(a,2),length(b))
        eid = 'Size:notEqual';
        msg = 'Second dimension of first input must equal length of second input.';
        throwAsCaller(MException(eid,msg))
    end
end

function mustBePercent(a)
   if  max(a < 0) && max(a > 100)
        eid = 'Val:notRange';
        msg = 'Noise must be in percent range o to 100';
        throwAsCaller(MException(eid,msg))
    end
end
