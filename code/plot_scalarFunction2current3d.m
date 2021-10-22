function plot_scalarFunction2current3d(x,y,z, valFunction, fileName, plot_program)
    arguments
        x (:,:) double {mustBeReal}
        y (:,:) double {mustBeReal}
        z (:,:) double {mustBeReal}
        valFunction (3,:,:) double {mustBeReal, mustBeDimensionFit(valFunction,x,y,z)}
        fileName (1,:) char = 'function'
        plot_program (1,:) char {mustBeMember(plot_program, {'none','matlab', 'tikz'})} = 'tikz'
    end
    
    fileNameShort = fileName(find(fileName == filesep,1,'last')+1:find(fileName == '.',1,'last')-1);
    plotFile = struct('type', fileNameShort, 'x', x, 'y', y, 'z', z, 'valSolution', valFunction);    
    switch plot_program
        case 'matlab'
            figure;
            title_name = replace(fileNameShort,'_',' ');
            valNorm = squeeze(sqrt(sum(valFunction.^2,1)));
            surf(x, y, z, valNorm, 'EdgeColor', 'none');
            hold on;
            quiver3(x, y, z, squeeze(valFunction(1,:,:))./valNorm, squeeze(valFunction(2,:,:))./valNorm, squeeze(valFunction(3,:,:))./valNorm);
            shading interp; colorbar; title(title_name); axis equal;
            view(-65,15);
        case 'tikz'
            fileNameSuffix = fileName(1:find(fileName == '.',1,'last')-1);
            fileNameQuiver = [fileNameSuffix, '_quiver', '.tex'];
    
            evalRadius = max(sqrt(x(:).^2+y(:).^2+z(:).^2));
            valNorm = squeeze(sqrt(sum(valFunction.^2,1)));
            
            plotExportPGFSphereQuiver(fileNameQuiver, 'Bilder/',  [repmat([-evalRadius, evalRadius], 3,1); 0, max(max(squeeze(valNorm)))], [{'$x$'},{'$y$'},{'$z$'}], [-65,15], true, x, y, z, valFunction);
        case 'none'
            disp('Plotting is disabled.');
    end
      
    save(['../solution/', fileNameShort, '_quiver.mat'], '-struct', 'plotFile');   
end

function mustBeDimensionFit(v,x,y,z)
    if ~isequal(size(x),size(y)) || ~isequal(size(z),size(y)) || ~isequal(size(x),size(squeeze(v(1,:,:))))
        eid = 'Size:notEqual';
        msg = 'Second dimension of first input must equal length of second input.';
        throwAsCaller(MException(eid,msg))
    end
end    