function plotExportPGFSphereLog(xData, yData, zData, cData, fileName, viewAngles, axisLimits, labels, hasColorbar, tikzSubPath)
    % Original function written by Samuel Leweke in 2017
    % Improved and adapted in 2021 by Sarah Leweke
    arguments
       xData (:,:) double {mustBeReal}
       yData (:,:) double {mustBeReal}
       zData (:,:) double {mustBeReal}
       cData (:,:) double {mustBeReal} = zData
       fileName (1,:) char = 'plotOnSphere.tex'
       viewAngles (1,2) double {mustBeReal} = [-65, 15]
       axisLimits (:,2) double {mustBeReal, mustBe34Dim(axisLimits)} = [-1, 1; -1, 1;-1,1]*0.074
       labels (1,:) cell {mustBe3Dim(labels)} = {'x','y','z'}
       hasColorbar (1,1) logical = true
       tikzSubPath (1,:) char = 'Bilder/'
    end
    
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
    
    % Downsample via interpolation if necessary; more than 60*60 points are
    % not compilable wiht pdflatex (2018) // out of memory
    nPoints = 60;

    % Threshold factor for detecting negative color data values (0.1% of maximum positive value)
    negativeColorThreshold = 0.001;

    if numel(xData) > nPoints * nPoints

        % Generate new grid
        phi = linspace(0,2*pi,nPoints);
        theta = linspace(0,pi/2,nPoints);
        [PhiNew, ThetaNew] = meshgrid(phi, theta);

        newX = radius.*sqrt(1-cos(ThetaNew).^2).*cos(PhiNew);
        newY = radius.*sqrt(1-cos(ThetaNew).^2).*sin(PhiNew);
        newZ = radius.*cos(ThetaNew);

        if min(min(cData == zData))
            F = scatteredInterpolant(xData(:), yData(:), zData(:), 'natural');
            newC = F(newX, newY);
            newZ = newC;
        else
            F = scatteredInterpolant(xData(:), yData(:), zData(:), cData(:), 'nearest');
            newC = F(newX, newY, newZ);
        end

        xData = newX;
        yData = newY;
        zData = newZ;
        cData = newC;
    end

    %  Export data
    if min(min(cData == zData))
        data = [ xData(:) yData(:) zData(:) ];
        % 8 digits, '-ascii','-double' for 16 digits 
        save([fileNameWithoutExt, '.dat'], 'data', '-ascii');
        pointMetaSrc = '';
        metaIndex = '';
%        zLim = cLim;
    else
        data = [ xData(:) yData(:) zData(:) cData(:) ];
        % 8 digits, '-ascii','-double' for 16 digits 
        save([fileNameWithoutExt, '.dat'], 'data', '-ascii');
        pointMetaSrc = 'point meta=explicit,';
        metaIndex = 'meta expr=ln(\thisrowno{3}),';

        if (size(axisLimits, 1) < 4)
            % Only overwrite color axis if limits are not given explicitly
            axisLimits(4,:) = [min(cData(:)), max(cData(:))];

            % Clip lower color limit to 0 if not sufficiently negative
%             if axisLimits(4,1) > -negativeColorThreshold * axisLimits(4,2)
%                 axisLimits(4,1) = 0;
%             end
        end
    end
    
    fid = fopen(fileName, 'w');
    fprintf(fid, '\\tikzsetnextfilename{%s}\n', baseName);
    fprintf(fid, '\\begin{tikzpicture}\n%%\n\\begin{axis}[%%\n');
    fprintf(fid, '\tstandardStyleSphere,\n\twidth=\\figurewidth,\n\theight=\\figurewidth,\n');
    fprintf(fid, '\tview={%g}{%g},\n', viewAngles);
    fprintf(fid, '\txmin=%g,\n\txmax=%g,\n', axisLimits(1,:));
    fprintf(fid, '\tymin=%g,\n\tymax=%g,\n', axisLimits(2,:));
    fprintf(fid, '\tzmin=%g,\n\tzmax=%g,\n', axisLimits(3,:));
    fprintf(fid, '\txlabel={%s},\n', labels{1});
    fprintf(fid, '\tylabel={%s},\n', labels{2});
    fprintf(fid, '\tzlabel={%s},\n', labels{3});
    if ~min(min(cData == zData)) && (size(axisLimits, 1) >= 4)
        fprintf(fid, '\tpoint meta min=%g,\n\tpoint meta max=%g,\n', [floor(log10(axisLimits(4,1))), ceil(log10(axisLimits(4,2)))]);
    end
    if hasColorbar

        if (~min(min(cData == zData)))
            % Check for negative values
            if axisLimits(4,1) <= -negativeColorThreshold * axisLimits(4,2)
                % We have some serious negative values

                % Calculate starting point of colormap
                lowerPoint = floor(-630 * (-axisLimits(4,1)) / axisLimits(4,2));
                % Emit custom colormap based on parula
%                fprintf(fid, '\tcolormap={changepoint}{[1pt]\nrgb255(%dpt)=(28,23,69)\nrgb255(-1pt)=(25,25,122)\n', lowerPoint);
                fprintf(fid, '\tcolormap={changepoint}{[1pt]\nrgb255(%dpt)=(28,23,69)\n', lowerPoint);
                fprintf(fid, 'rgb(0pt)=(0.2081,0.1663,0.5292)\nrgb(10pt)=(0.2116,0.1898,0.5777)\nrgb(20pt)=(0.2123,0.2138,0.627)\nrgb(30pt)=(0.2081,0.2386,0.6771)\nrgb(40pt)=(0.1959,0.2645,0.7279)\nrgb(50pt)=(0.1707,0.2919,0.7792)\nrgb(60pt)=(0.1253,0.3242,0.8303)\n');
                fprintf(fid, 'rgb(70pt)=(0.0591,0.3598,0.8683)\nrgb(80pt)=(0.0117,0.3875,0.882)\nrgb(90pt)=(0.006,0.4086,0.8828)\nrgb(100pt)=(0.0165,0.4266,0.8786)\nrgb(110pt)=(0.0329,0.443,0.872)\nrgb(120pt)=(0.0498,0.4586,0.8641)\nrgb(130pt)=(0.0629,0.4737,0.8554)\n');
                fprintf(fid, 'rgb(140pt)=(0.0723,0.4887,0.8467)\nrgb(150pt)=(0.0779,0.504,0.8384)\nrgb(160pt)=(0.0793,0.52,0.8312)\nrgb(170pt)=(0.0749,0.5375,0.8263)\nrgb(180pt)=(0.0641,0.557,0.824)\nrgb(190pt)=(0.0488,0.5772,0.8228)\nrgb(200pt)=(0.0343,0.5966,0.8199)\n');
                fprintf(fid, 'rgb(210pt)=(0.0265,0.6137,0.8135)\nrgb(220pt)=(0.0239,0.6287,0.8038)\nrgb(230pt)=(0.0231,0.6418,0.7913)\nrgb(240pt)=(0.0228,0.6535,0.7768)\nrgb(250pt)=(0.0267,0.6642,0.7607)\nrgb(260pt)=(0.0384,0.6743,0.7436)\nrgb(270pt)=(0.059,0.6838,0.7254)\n');
                fprintf(fid, 'rgb(280pt)=(0.0843,0.6928,0.7062)\nrgb(290pt)=(0.1133,0.7015,0.6859)\nrgb(300pt)=(0.1453,0.7098,0.6646)\nrgb(310pt)=(0.1801,0.7177,0.6424)\nrgb(320pt)=(0.2178,0.725,0.6193)\nrgb(330pt)=(0.2586,0.7317,0.5954)\nrgb(340pt)=(0.3022,0.7376,0.5712)\n');
                fprintf(fid, 'rgb(350pt)=(0.3482,0.7424,0.5473)\nrgb(360pt)=(0.3953,0.7459,0.5244)\nrgb(370pt)=(0.442,0.7481,0.5033)\nrgb(380pt)=(0.4871,0.7491,0.484)\nrgb(390pt)=(0.53,0.7491,0.4661)\nrgb(400pt)=(0.5709,0.7485,0.4494)\nrgb(410pt)=(0.6099,0.7473,0.4337)\n');
                fprintf(fid, 'rgb(420pt)=(0.6473,0.7456,0.4188)\nrgb(430pt)=(0.6834,0.7435,0.4044)\nrgb(440pt)=(0.7184,0.7411,0.3905)\nrgb(450pt)=(0.7525,0.7384,0.3768)\nrgb(460pt)=(0.7858,0.7356,0.3633)\nrgb(470pt)=(0.8185,0.7327,0.3498)\nrgb(480pt)=(0.8507,0.7299,0.336)\n');
                fprintf(fid, 'rgb(490pt)=(0.8824,0.7274,0.3217)\nrgb(500pt)=(0.9139,0.7258,0.3063)\nrgb(510pt)=(0.945,0.7261,0.2886)\nrgb(520pt)=(0.9739,0.7314,0.2666)\nrgb(530pt)=(0.9938,0.7455,0.2403)\nrgb(540pt)=(0.999,0.7653,0.2164)\nrgb(550pt)=(0.9955,0.7861,0.1967)\n');
                fprintf(fid, 'rgb(560pt)=(0.988,0.8066,0.1794)\nrgb(570pt)=(0.9789,0.8271,0.1633)\nrgb(580pt)=(0.9697,0.8481,0.1475)\nrgb(590pt)=(0.9626,0.8705,0.1309)\nrgb(600pt)=(0.9589,0.8949,0.1132)\nrgb(610pt)=(0.9598,0.9218,0.0948)\nrgb(620pt)=(0.9661,0.9514,0.0755)\n');
                fprintf(fid, 'rgb(630pt)=(0.9763,0.9831,0.0538)\n},\n');
            end
        end

        fprintf(fid, '\tcolorbar=true,\n');
        fprintf(fid, '\tcolorbar horizontal,\n'); 
        fprintf(fid, '\tcolorbar style={at={(0.5,1.03)},anchor=south, xticklabel pos=upper, xticklabel=$10^{\\pgfmathprintnumber{\\tick}}$},\n');
    end
    fprintf(fid, ']\n');
    fprintf(fid, '\\addplot3[surf,shader=interp,z buffer=sort,row sep=crcr,mesh/ordering=rowwise,mesh/cols=%d,%s] table[x index=0,y index=1,z index=2,%sheader=false] {%s%s.dat};\n', size(xData, 2), pointMetaSrc, metaIndex, tikzSubPath, baseName);
    fprintf(fid, '\n\\end{axis}\n\\end{tikzpicture}%%');
    fclose(fid);
end

function mustBe34Dim(a)
    if size(a,1) < 3 && size(a,1) > 4
        eid = 'Size:notPermitted';
        msg = 'Axis Dimension for three dimensional plot must be 3 or 4 (data) dimensions';
        throwAsCaller(MException(eid,msg))
    end
end

function mustBe3Dim(a)
    if length(a) ~= 3
        eid = 'Size:notPermitted';
        msg = 'Axis Dimension must be of length three (3d plot)';
        throwAsCaller(MException(eid,msg))
    end
end
