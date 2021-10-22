function plot_scalarFunction3d(x,y,z, valFunction, fileName, plot_program, doLog)
    arguments
        x (:,:) double {mustBeReal}
        y (:,:) double {mustBeReal}
        z (:,:) double {mustBeReal}
        valFunction (:,:) double {mustBeReal, mustBeDimensionFit(valFunction,x,y,z)}
        fileName (1,:) char = 'function'
        plot_program (1,:) char {mustBeMember(plot_program, {'none','matlab', 'tikz'})} = 'tikz'
        doLog (1,1) logical = false
    end
    
    fileNameShort = fileName(find(fileName == filesep,1,'last')+1:find(fileName == '.',1,'last')-1);
    plotFile = struct('type', fileNameShort, 'x', x, 'y', y, 'z', z, 'valSolution', valFunction);    
    switch plot_program
        case 'matlab'
            figure;
            title_name = replace(fileNameShort,'_',' ');
            surf(x,y,z,valFunction); shading interp; colorbar; title(title_name)
        case 'tikz'
            if doLog
                plotExportPGFSphereLog(x, y, z, valFunction, fileName)
            else
                plotExportPGFSphere(x, y, z, valFunction, fileName)
            end
        case 'none'
            disp('Plotting is disabled.')
    end
      
    save(['../solution/', fileNameShort, '.mat'], '-struct', 'plotFile');   
end

function mustBeDimensionFit(v,x,y,z)
    if ~isequal(size(x),size(y)) || ~isequal(size(z),size(y)) || ~isequal(size(x),size(v))
        eid = 'Size:notEqual';
        msg = 'Second dimension of first input must equal length of second input.';
        throwAsCaller(MException(eid,msg))
    end
end    