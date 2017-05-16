function [pVsDist,dBins,frameTimes]  = migrationParamVsDistance(tracks,tPos,frameTimes,fieldName,varargin)

%Makes map of average value across all tracks of specified migration
%parameter vs distance and time
%
% Input: 
%   tracks - structure with tracks
%   tPos - target position
%   frameTime - time corresponding to each frame.
%   fieldName - name of field in track structure to display
%



%% ---------- Input ---------- %%

ip = inputParser;
ip.addParameter('BinWidth',10,@(x)(isscalar(x) && x > 0));%Distance bin size (in same units as track positions, usually microns)
ip.addParameter('MaxDistance',5e2,@(x)(isscalar(x) && x > 0));%Maximum distance to include (in same units as track positions, usually microns)
ip.addParameter('MinDistance',0,@(x)(isscalar(x) && x >= 0));%Minimum distance to include (in same units as track positions, usually microns)
ip.addParameter('DisplayName',fieldName,@ischar);%Display name for field to use in figure. Default is same as fieldName.
ip.addParameter('ShowFigures',true,@islogical);%If true, display figures
ip.addParameter('FillNaNs',true,@islogical);%If true, interpolation will be used to fill NaN areas (e.g. where no cells were present). Should only be a small number of bins, if more than 1% error is thrown.
ip.addParameter('BatchMode',false,@islogical);%If true, figure display is supressed for batch processing.
ip.parse(varargin{:});
p = ip.Results;

%% ---------- Init ---------- %%


dBins = p.MinDistance:p.BinWidth:p.MaxDistance;
nBins = numel(dBins)-1;
nFrames = numel(frameTimes);

iFirst = arrayfun(@(x)(min(x.Frame)),tracks);
iLast = arrayfun(@(x)(max(x.Frame)),tracks);


nPerDist = nan(nFrames,nBins);
pVsDist = nan(nFrames,nBins);



%% ------ Make Binned Averaged Map --------- %%

%TEMP - vectorize this shit?? kd-tree it? fast enough for now...
for j = 1:nFrames
    
    %Get tracks which exist in current frame
    iExist = find(iFirst <= j & iLast >= j);
    nCurr = numel(iExist);
    pos = nan(nCurr,2);
    par = nan(nCurr,2);
    for k = 1:nCurr
        iFr = find(tracks(iExist(k)).Frame == j);
        if ~isempty(iFr) %check for gaps...
            %Get the position and the parameter value at this frame
            pos(k,:) = [tracks(iExist(k)).x(iFr) tracks(iExist(k)).y(iFr)];                
            if ~isempty(tracks(iExist(k)).(fieldName)) && numel(tracks(iExist(k)).(fieldName))>= iFr%Parameters may not be defined for every track / timepoint
                par(k) = tracks(iExist(k)).(fieldName)(iFr);
            end
        end
    end
    
    %Get distance to target
    d = sqrt(sum( bsxfun(@minus,pos,tPos) .^2,2));
    %Get number vs distance histogram
    [nPerDist(j,:),~,binInd] = histcounts(d,dBins);
    
    
    %Get avg param vs distance 
    pVsDist(j,:) = accumarray(binInd(binInd>0),par(binInd>0),[nBins 1],@nanmean);
    
    
end

if p.FillNaNs && nnz(isnan(pVsDist)) > 0
    %Make sure we don't have too many missing values as this likely
    %indicates a problem with the data, bin size etc.
    misVals = isnan(pVsDist);
    if (nnz(misVals) / numel(pVsDist)) > .01
        error('Too many NaN values to interpolate! Aborting!')
    end
    %Fill in these values with linear interpolation
    [X,Y] = meshgrid(1:nBins,1:nFrames);
    intFun = scatteredInterpolant(X(~misVals),Y(~misVals),pVsDist(~misVals));    
    pVsDist(misVals) = intFun(X(misVals),Y(misVals));        
    
end
    

pVsDist = pVsDist';%Transpose to match with functions for fitting and modelling. Yeah I know, iknow, it's just easier to fix it here..

if p.ShowFigures
    
    if ~p.BatchMode        
        cf = figure;
    else
        cf = figure('Visible','off');
    end
    imHan = imagesc(frameTimes,dBins,pVsDist);
    ca = get(cf,'CurrentAxes');
    set(ca,'YDir','normal')
    xlabel('Time, seconds')
    ylabel('Distance, microns')
    cHan = colorbar;
    cHan.Label.String = p.DisplayName;
    
end

