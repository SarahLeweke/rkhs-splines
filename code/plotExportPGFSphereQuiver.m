function plotExportPGFSphereQuiver(fileName, tikzSubPath, axisLimits, labels, viewAngles, hasColorbar, xData, yData, zData, valData, disableDownsampling)
%PLOTEXPORTPGFSPHEREQUIVER Exports a sphere + quiver plot (on top) to pgfplots
%
% plotExportPGFSphereQuiver('bla2.tex', '', [repmat([-0.6, 0.6], 3,1); 0, 1], [{'X'},{'Y'},{'Z'}], [-130,26], true, evalPos.evalRadius*evalPos.evalPos1, evalPos.evalRadius*evalPos.evalPos2, evalPos.evalRadius*evalPos.evalPos3, val);

	if (length(size(valData)) ~= 3) || (size(valData, 1) ~= 3)
        error('valData has to be of size 3 x nGrid x nGrid');
	end

	if nargin <= 10
		disableDownsampling = false;
	end
	
    % Use only every nQuiverPoints-point for quiver arrows (default: every
    % 7th point)
    nQuiverPoints = 5;

    % Number of downsampling points
    nPoints = 61;
    
    % Hidden arrow (backside) removal tolerance
    tolHiddenPoints = 1e-2;
    
    % Tolerance for detecting arrows that point inside the sphere
    tolOutsidePointingArrows = 1e-1;
    
    % Extract file name
    idx = find(fileName == filesep,1,'last');
    idxDot = find(fileName == '.',1,'last') - 1;
    if ~isempty(idx)
        baseName = fileName(idx+1:idxDot);
    else
        baseName = fileName(1:idxDot);
    end

    fileNameWithoutExt = fileName(1:idxDot);
    
    if ~isempty(tikzSubPath)
        if tikzSubPath(end) ~= filesep
            tikzSubPath = [tikzSubPath filesep];
        end
    end

    % Compute radius
    radius = sqrt(xData.^2 + yData.^2 + zData.^2);
    radius = max(radius(:));
    
    % Downsample via interpolation if necessary
    if (numel(xData) > 2 * nPoints * nPoints) && ~disableDownsampling

        % Generate new grid
        phi = linspace(0,2*pi,2*nPoints);
        theta = linspace(0,pi/2,nPoints);
        [PhiNew, ThetaNew] = meshgrid(phi, theta);

        newX = radius.*sqrt(1-cos(ThetaNew).^2).*cos(PhiNew);
        newY = radius.*sqrt(1-cos(ThetaNew).^2).*sin(PhiNew);
        newZ = radius.*cos(ThetaNew);
        
        newVal = zeros(3, size(newX,1), size(newX,2));
        
        % Interpolate each component separately
        for i = 1:3
%                F = scatteredInterpolant(xData(:), yData(:), zData(:), cData(:), 'natural');
            curVal = squeeze(valData(i,:,:));
            F = scatteredInterpolant(xData(:), yData(:), zData(:), curVal(:), 'nearest');
            newC = F(newX, newY, newZ);
%                F = griddedInterpolant(ThetaOld, PhiOld, cData, 'spline');
%                newC = F(ThetaNew, PhiNew);
            newVal(i,:,:) = newC;
        end
        xData = newX;
        yData = newY;
        zData = newZ;
        valData = newVal;
    end

    meshCols = size(xData, 2);

    %  Export data for sphere
    cValData = squeeze(sqrt(sum(valData.^2,1)));
    data = [ xData(:) yData(:) zData(:) cValData(:)];
    % 8 digits, '-ascii','-double' for 16 digits
    save([fileNameWithoutExt, '-sphere.dat'], 'data', '-ascii');

    % Export data for quiver
    xValData = squeeze(valData(1,:,:));
    xValData = xValData(:);
    yValData = squeeze(valData(2,:,:));
    yValData = yValData(:);
    zValData = squeeze(valData(3,:,:));
    zValData = zValData(:);
    xData = xData(:);
    yData = yData(:);
    zData = zData(:);
    data = [ xData(1:nQuiverPoints:end) yData(1:nQuiverPoints:end) zData(1:nQuiverPoints:end) xValData(1:nQuiverPoints:end) yValData(1:nQuiverPoints:end) zValData(1:nQuiverPoints:end)];

    % Compute normal vector facing towards us
    T = viewmtx(viewAngles(1), viewAngles(2));
    T = T(1:3, 1:3); % Remove useless 4th component
    normVec = [0,0,1] * T; % normVec has norm = 1

    normArrowBase = sqrt(sum(data(:,1:3).^2, 2));

    % Keep only points that are on the frontside of the sphere (from our
    % point of view)
    idxFront = data(:,1:3) * normVec' >= -tolHiddenPoints * normArrowBase;

    % Keep only arrows that point outside the sphere
    idxDir = diag(data(:,4:6) * data(:,1:3)') >= -tolOutsidePointingArrows .* sqrt(sum(data(:,4:6).^2, 2)) .* normArrowBase;

    arrowData = data(idxFront & idxDir, :);

    % Calculate scaling factor for arrows
    arrowData(:,4:end) = arrowData(:,4:end) .* autoScaleArrows(arrowData(:,1), arrowData(:,2), arrowData(:,3), arrowData(:,4:end));
    arrowLen = sqrt(sum(arrowData(:,4:end).^2, 2));
    if abs(max(arrowLen)) < eps
        arrowData = [arrowData, arrowLen];
    else
        arrowData = [arrowData, arrowLen ./ max(arrowLen)]; % Scale arrow Length to [0, 1] range
    end
    
    % 8 digits, '-ascii','-double' for 16 digits
    save([fileNameWithoutExt, '-quiver.dat'], 'arrowData', '-ascii');

    % Export points on front that point inwards
    %   When arrows point inwards, we cannot plot the arrow. Instead we
    %   plot the tail of the arrow as point on the surface
    pointData = data(idxFront & ~idxDir, 1:3);
    if ~isempty(pointData)
        % 8 digits, '-ascii','-double' for 16 digits
        save([fileNameWithoutExt, '-points.dat'], 'pointData', '-ascii');
    end
    
    fid = fopen(fileName, 'w');
    fprintf(fid, '\\tikzsetnextfilename{%s}\n', baseName);
    fprintf(fid, '\\begin{tikzpicture}\n%%\n\\begin{axis}[%%\n');
    fprintf(fid, '\tstandardStyleSphereQuiver,\n\twidth=\\figurewidth,\n\theight=\\figurewidth,\n');
    fprintf(fid, '\tview={%g}{%g},\n', viewAngles);
  fprintf(fid, '\txmin=%g,\n\txmax=%g,\n', 100*axisLimits(1,:));
    fprintf(fid, '\tymin=%g,\n\tymax=%g,\n', 100*axisLimits(2,:));
    fprintf(fid, '\tzmin=%g,\n\tzmax=%g,\n', 100*axisLimits(3,:));
    fprintf(fid, '\txlabel={$%s$},\n', labels{1});
    fprintf(fid, '\tylabel={$%s$},\n', labels{2});
    fprintf(fid, '\tzlabel={$%s$},\n', labels{3});
    fprintf(fid, '\tx unit={\\si{\\centi\\metre}},\n');
    fprintf(fid, '\ty unit={\\si{\\centi\\metre}},\n');
    fprintf(fid, '\tz unit={\\si{\\centi\\metre}},\n');
	fprintf(fid, '\txlabel style={font=\\footnotesize,sloped},\n');
	fprintf(fid, '\tylabel style={font=\\footnotesize,sloped},\n');
	fprintf(fid, '\tzlabel style={font=\\footnotesize,sloped},\n');
    if (size(axisLimits, 1) >= 4)
        fprintf(fid, '\tpoint meta min=%g,\n\tpoint meta max=%g,\n', axisLimits(4,:));
    end
    if hasColorbar
        fprintf(fid, '\tcolorbar=true,\n');
        fprintf(fid, '\tcolorbar horizontal,\n'); 
        fprintf(fid, '\tcolorbar style={at={(0.5,1.03)},anchor=south, xticklabel pos=upper,}\n');
    end
    fprintf(fid, ']\n');
    fprintf(fid, '\\addplot3[surf,shader=interp,z buffer=sort,mesh/ordering=rowwise,mesh/cols=%d,point meta=explicit,] table[x expr={100*\\thisrowno{0}},y expr={100*\\thisrowno{1}},z expr={100*\\thisrowno{2}},meta index=3,header=false] {%s%s-sphere.dat};\n', meshCols, tikzSubPath, baseName);
    if ~isempty(pointData)
        fprintf(fid, '\\addplot3[black,only marks,every mark/.append style={scale=0.1}] table[x expr={100*\\thisrowno{0}},y expr={100*\\thisrowno{1}},z expr={100*\\thisrowno{2}},header=false] {%s%s-points.dat};\n', tikzSubPath, baseName);
    end
    fprintf(fid, '\\addplot3[black,visualization depends on=\\thisrowno{6}\\as\\localmeta,quiver={u={100*\\thisrowno{3}},v={100*\\thisrowno{4}},w={100*\\thisrowno{5}},scale arrows=1.0,every arrow/.append style={-{Stealth[round,scale={0.75*max(0.001,\\localmeta)}]},},},] table[x expr={100*\\thisrowno{0}},y expr={100*\\thisrowno{1}},z expr={100*\\thisrowno{2}},header=false] {%s%s-quiver.dat};\n', tikzSubPath, baseName);
    fprintf(fid, '\n\\end{axis}\n\\end{tikzpicture}%%');
    fclose(fid);
end

function [autoscaleFactor] = autoScaleArrows(xData, yData, zData, valData)
    % Base autoscale value on average spacing in the x and y
    % directions.  Estimate number of points in each direction as
    % either the size of the input arrays or the effective square
    % spacing if x and y are vectors.
    n=sqrt(numel(xData));
    m=n;
    delx = diff([min(xData(:)) max(xData(:))])/n;
    dely = diff([min(yData(:)) max(yData(:))])/m;
    delz = diff([min(zData(:)) max(zData(:))])/max(m,n);
    del = sqrt(delx.^2 + dely.^2 + delz.^2);
    if del>0
        len = squeeze(sqrt(sum(valData.^2,2))) ./ del;
        maxlen = max(len(:));
    else
        maxlen = 0;
    end
    
    if maxlen>0
        autoscaleFactor = 0.9 / maxlen;
    else
        autoscaleFactor = 0.9;
    end
end
