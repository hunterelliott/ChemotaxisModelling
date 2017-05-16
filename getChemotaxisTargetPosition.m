function [tPos,minD,maxD] = getChemotaxisTargetPosition(imFile,trk,iChan,figOutDir,batchMode)


if nargin < 3
    figOutDir = '';
end

%Create MD to simplify image/metadata reading
MD = bfImport(imFile);

% ----  Segment target and get centroid ---- %

%We assume stationary target and just use first frame
im = MD.channels_(iChan).loadImage(1,1);

imThresh = thresholdOtsu(im);

tMask = filterGauss2D(im,1.2) > imThresh;

tMask = bwLargestObj(tMask);

rp = regionprops(tMask,'Centroid','Area');
tPos = rp.Centroid * MD.pixelSize_/1e3;%Convert to microns

%Check image origin to adjust for stage position (imaris uses absolute not
%relative positions, taking stage position in metadate into account). For
%whatever reason bio-formats can't read this info so we assume tracks touch
%every boundary (always the case in our data).
imOrig = [min(vertcat(trk.x)) min(vertcat(trk.y))];
imOrig = imOrig - MD.pixelSize_/2e3; %Imaris uses pixel centers for track positions

%Adjus target position based on origin

tPos = tPos + imOrig;

% ---- Valid distance determination --- %
%We set a range of valid distances since images are square and we are
%approximating the target as a point

%Get closest distance to image border, again assuming tracks touch all
%boundaries
alldX = [(tPos(1) - vertcat(trk.x)) , (tPos(2) -vertcat(trk.y)) ];
minXtoEdge = min(abs(min(alldX(:,1))),max(alldX(:,1)));
minYtoEdge = min(abs(min(alldX(:,2))),max(alldX(:,2)));
maxD = min(minXtoEdge,minYtoEdge);


%Because we use point approximation, close distances will have significant
%error due to finite radius of target. So we specify a minimum distance
%where this error is acceptable.
maxDistErr = .75;%Maximum fractional distance error;
rTarget = sqrt(rp.Area/pi);%Equivalent target radius;
minD = 2*rTarget / maxDistErr;%Resulting minimum valid distance
minD = minD * MD.pixelSize_ / 1e3;


if ~isempty(figOutDir)
    if batchMode
        cf = figure('Visible','off');
    else        
        cf = figure;
    end
    imshow(im,[])
    %saturateImageColormap(cf.CurrentAxes);
    hold on
    cellfun(@(x)(plot(x(:,2),x(:,1),'r-')),bwboundaries(tMask))
    tPosPlot = tPos - imOrig;
    tPosPlot = tPosPlot / MD.pixelSize_  *1e3;
    plot(tPosPlot(1),tPosPlot(2),'rx')
    pR = minD / MD.pixelSize_ *1e3;
    rectangle('Position',[(tPosPlot -pR) ,pR*2 pR*2],'Curvature',[1 1],'FaceColor','none','EdgeColor','g','LineWidth',2);
    pR = maxD / MD.pixelSize_ *1e3;
    rectangle('Position',[(tPosPlot -pR) ,pR*2 pR*2],'Curvature',[1 1],'FaceColor','none','EdgeColor','g','LineWidth',2);
    
    mfFigureExport(cf,[figOutDir filesep 'target position'])
end





