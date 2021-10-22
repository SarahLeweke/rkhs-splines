function plotExportPGFSphereScatter(fileName, tikzSubPath, sensorPos, measValues, hFig)
    
    scaleMeasValuesToColorBar = false;

    % Extract file name
    idx = find(fileName == filesep,1,'last');
    idxDot = find(fileName == '.',1,'last') - 1;
    if ~isempty(idx)
        baseName = fileName(idx+1:idxDot);
    else
        baseName = fileName(1:idxDot);
    end

    fileName = fileName;
    fileNameWithoutExt = fileName(1:idxDot);
    
    if ~isempty(tikzSubPath)
        if tikzSubPath(end) ~= filesep
            tikzSubPath = [tikzSubPath filesep];
        end
    end
    
    if size(sensorPos, 2) ~= 3
        sensorPos = sensorPos';
    end
    
    % Extract data from figure
    if (nargin < 5) || isempty(hFig)
        hFig = gca;
    end
	axesObjs = get(hFig, 'Children');
    for i = 1:length(axesObjs)
        curObj = axesObjs(i);
        objType = get(curObj, 'Type');
        if strcmpi(objType, 'axes')
            % Axes: Inspect children
            dataObjs = get(curObj, 'Children');
            objType = get(dataObjs, 'Type');
            xData = get(dataObjs, 'XData');
            yData = get(dataObjs, 'YData');
            zData = get(dataObjs, 'ZData');
            cData = get(dataObjs, 'CData');
            
            % Transpose data
            xData = xData';
            yData = yData';
            zData = zData';
            cData = cData';

            % Get axis limits
            xLim = get(curObj, 'XLim');
            yLim = get(curObj, 'YLim');
            zLim = get(curObj, 'ZLim');
            cLim = get(curObj, 'CLim');
            viewAngles = get(curObj, 'View');
            
            labels = cellfun(@(x) get(x, 'String'), get(curObj,{'XLabel' 'YLabel' 'ZLabel'}), 'UniformOutput', false);
                        
            % Compute radius
            radius = sqrt(xData.^2 + yData.^2 + zData.^2);
            radius = max(radius(:));

            % Downsample via interpolation if necessary
            nPoints = 60;

            if numel(xData) > nPoints * nPoints

                % Generate old grid
                phi = linspace(0,2*pi,size(zData,1));
                theta = linspace(0,pi/2,size(zData,2));
                [PhiOld, ThetaOld] = meshgrid(phi, theta);

                % Generate new grid
                phi = linspace(0,2*pi,nPoints);
                theta = linspace(0,pi/2,nPoints);
                [PhiNew, ThetaNew] = meshgrid(phi, theta);

                newX = radius.*sqrt(1-cos(ThetaNew).^2).*cos(PhiNew);
                newY = radius.*sqrt(1-cos(ThetaNew).^2).*sin(PhiNew);
                newZ = radius.*cos(ThetaNew);

%                F = scatteredInterpolant(xData(:), yData(:), zData(:), cData(:), 'natural');
                F = scatteredInterpolant(xData(:), yData(:), zData(:), cData(:), 'nearest');
                newC = F(newX, newY, newZ);
%                F = griddedInterpolant(ThetaOld, PhiOld, cData, 'spline');
%                newC = F(ThetaNew, PhiNew);

                xData = newX;
                yData = newY;
                zData = newZ;
                cData = newC;
            end

            % Scale measurements to colorbar if required
            if scaleMeasValuesToColorBar
                % Translate to positive range
%                measValues = min(measValues) + measValues;
                % Compress to [0, cMax]
                measValues = measValues ./ max(measValues) * max(cLim);
            end
            
            % Update axis limits due to additional scatter points
            scatterXLim = [floor(min(sensorPos(:,1)) .* 10) ./ 10, ceil(max(sensorPos(:,1)) .* 10) ./ 10];
            scatterYLim = [floor(min(sensorPos(:,2)) .* 10) ./ 10, ceil(max(sensorPos(:,2)) .* 10) ./ 10];
            scatterZLim = [floor(min(sensorPos(:,3)) .* 10) ./ 10, ceil(max(sensorPos(:,3)) .* 10) ./ 10];
            xLim = [min(scatterXLim(1), xLim(1)), max(scatterXLim(2), xLim(2))];
            yLim = [min(scatterYLim(1), yLim(1)), max(scatterYLim(2), yLim(2))];
            zLim = [min(scatterZLim(1), zLim(1)), max(scatterZLim(2), zLim(2))];
            
            % Remove hidden points (i.e., point behind the sphere from
            % point of view)
            camPos = get(curObj, 'CameraPosition');
            device = load_device;
            radius = device.radiusCerebrum;
            camPos = [1.5054,1.535,0.5854];%[-1.6657,-1.0606,0.5854];
            
            idxRemove = false(size(sensorPos, 1), 1);
            for j = 1:size(sensorPos, 1)
                idxRemove(j) = testPointBehindSphere(camPos, radius, sensorPos(j,:));
            end
            measValues = measValues(:);
            data = [ sensorPos(~idxRemove, :) measValues(~idxRemove) ];
            save(fileName, 'data', '-ascii');
            
            
            measLim = [floor(min(measValues)), ceil(max(measValues))];
            
            %  Export data for sphere
            data = [ xData(:) yData(:) zData(:) cData(:) ];
            % 8 digits, '-ascii','-double' for 16 digits 
            save([fileNameWithoutExt, '-sphere.dat'], 'data', '-ascii');

            %  Export data for sensor points
            data = [ sensorPos(~idxRemove, :) measValues(~idxRemove) ];
            % 8 digits, '-ascii','-double' for 16 digits 
            save([fileNameWithoutExt, '-sensor.dat'], 'data', '-ascii');
            
            % Write tikz file
            fid = fopen(fileName, 'w');
            fprintf(fid, '\\tikzsetnextfilename{%s}\n', baseName);
            fprintf(fid, '\\begin{tikzpicture}\n');
            fprintf(fid, '\\pgfplotsset{colorbar/draw/.code={%%\n\t\\axis[every colorbar,colorbar shift,colorbar=false,at={(0.5,1.07)},anchor=south,xticklabel pos=upper,axis x line*=top,]\n\t\t\\addplot graphics {};\n\t\\endaxis\n\t},\n}\n');
            fprintf(fid, '\\begin{axis}[%%\n');
            fprintf(fid, '\tstandardStyleSphere,\n\twidth=\\figurewidth,\n\theight=\\figurewidth,\n');
            fprintf(fid, '\tview={%g}{%g},\n', viewAngles);
            fprintf(fid, '\txmin=%g,\n\txmax=%g,\n', xLim);
            fprintf(fid, '\tymin=%g,\n\tymax=%g,\n', yLim);
            fprintf(fid, '\tzmin=%g,\n\tzmax=%g,\n', zLim);
            fprintf(fid, '\txlabel={%s},\n', labels{1});
            fprintf(fid, '\tylabel={%s},\n', labels{2});
            fprintf(fid, '\tzlabel={%s},\n', labels{3});
            fprintf(fid, '\tpoint meta min=%g,\n\tpoint meta max=%g,\n', cLim);
            fprintf(fid, '\tcolorbar=true,\n');
            fprintf(fid, '\tcolorbar horizontal,\n'); 
            fprintf(fid, ']\n');
            fprintf(fid, '\\addplot3[surf,shader=interp,z buffer=sort,mesh/ordering=rowwise,mesh/cols=%d,point meta=explicit,] table[x index=0,y index=1,z index=2,meta index=3,header=false] {%s%s-sphere.dat};\n', size(xData, 2), tikzSubPath, baseName);
            fprintf(fid, '\\end{axis}\n');
            fprintf(fid, '\\pgfplotsset{colorbar/draw/.code={%%\n\t\\axis[every colorbar,colorbar shift,colorbar=false,at={(0.5,1.07)},anchor=south,xticklabel pos=lower,axis x line*=bottom,axis y line=none,]\n\t\\endaxis\n\t},\n}\n');
            fprintf(fid, '\\begin{axis}[%%\n');
            fprintf(fid, '\tstandardStyleSphere,\n\twidth=\\figurewidth,\n\theight=\\figurewidth,\n');
            fprintf(fid, '\tview={%g}{%g},\n', viewAngles);
            fprintf(fid, '\txmin=%g,\n\txmax=%g,\n', xLim);
            fprintf(fid, '\tymin=%g,\n\tymax=%g,\n', yLim);
            fprintf(fid, '\tzmin=%g,\n\tzmax=%g,\n', zLim);
            fprintf(fid, '\tpoint meta min=%g,\n\tpoint meta max=%g,\n', measLim);
            fprintf(fid, '\taxis y line=none,\n\taxis x line=none,\n\taxis z line=none,\n');
            fprintf(fid, '\tcolorbar=true,\n');
            fprintf(fid, '\tcolorbar horizontal,\n'); 
            fprintf(fid, ']\n');
            fprintf(fid, '\\addplot3[z buffer=sort,scatter,scatter src=\\thisrowno{3},scatter/use mapped color={fill=mapped color,draw=mapped color},only marks,mark size=1pt,point meta=explicit] table[x index=0,y index=1,z index=2,meta index=3,header=false] {%s%s-sensor.dat};\n', tikzSubPath, baseName);
            fprintf(fid, '\n\\end{axis}\n');
            fprintf(fid, '\\pgfplotsset{colorbar/draw/.code={%%\n\t\\axis[every colorbar,colorbar shift,colorbar=false,]\n\t\t\\addplot graphics {};\n\t\\endaxis\n\t},\n}\n');
            fprintf(fid, '\\end{tikzpicture}%%');
            fclose(fid);
        end
    end
end

function retVal = testPointBehindSphere(viewPos, sphereRad, pt)
    
    retVal = false;

    % Obtain unit direction vector from camera to point
    dirVec = pt(:) - viewPos(:);
    dirVec = dirVec ./ norm(dirVec);
    dirVec = dirVec(:);

    % Calculate intersection position on ray
    % via https://en.wikipedia.org/wiki/Line?sphere_intersection
    temp = dirVec' * viewPos(:);
    detVal = temp^2 - sum(viewPos.^2) + sphereRad^2;
    
    % Check if line does not intersect sphere
    if detVal <= 0.0
        return;
    end
    
    % Compute first intersection (distance from view point)
    colPos = min(-temp + sqrt(detVal), -temp - sqrt(detVal));
    if colPos < norm(pt(:) - viewPos(:))
        % First intersection occurs earlier than distance to point
        % => Behind sphere
        retVal = true;
    end
    
end
