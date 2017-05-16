function [fracResPerDist,nResPerDist,nPerDist,dBins]  = migratingFractionVsDistance(tracks,iStart,tPos,frameTimes,varargin)

%% ---------- Input ---------- %%

ip = inputParser;
ip.addParameter('BinWidth',10,@(x)(isscalar(x) && x > 0));%Distance bin size (in same units as track positions, usually microns)
ip.addParameter('MaxDistance',300,@(x)(isscalar(x) && x > 0));%Maximum distance to include (in same units as track positions, usually microns)

ip.parse(varargin{:});
p = ip.Results;

%% ---------- Init ---------- %%


dBins = 0:p.BinWidth:p.MaxDistance;
nBins = numel(dBins);
nFrames = numel(frameTimes);

iFirst = arrayfun(@(x)(min(x.Frame)),tracks);
iLast = arrayfun(@(x)(max(x.Frame)),tracks);


nPerDist = nan(nFrames,nBins);
nResPerDist = nan(nFrames,nBins);

nPerDistCum = nan(nFrames,nBins);
nResPerDistCum = nan(nFrames,nBins);


%% ------ Calc --------- %%

%TEMP - vectorize this shit?? kd-tree it? fast enough for now...
for j = 1:nFrames
    
    %Get tracks which exist in current frame
    iExist = find(iFirst <= j & iLast >= j);
    nCurr = numel(iExist);
    pos = nan(nCurr,2);
    for k = 1:nCurr
        iFr = find(tracks(iExist(k)).Frame == j);
        if ~isempty(iFr) %check for gaps...
            pos(k,:) = [tracks(iExist(k)).x(iFr) tracks(iExist(k)).y(iFr)];                
        end
    end
    
    %Get distance to target
    d = sqrt(sum( bsxfun(@minus,pos,tPos) .^2,2));
    %Get number vs distance histogram
    nPerDist(j,:) = histc(d,dBins);
    nPerDistCum(j,:) = cumsum(nPerDist(j,:));
    
    %Get number responding vs distance histogram
    isRes = iStart(iExist) <= j;
    nResPerDist(j,:) = histc(d(isRes),dBins);
    nResPerDistCum(j,:) = cumsum(nResPerDist(j,:));
    
end

fracResPerDist = nResPerDist ./ nPerDist;

fracResPerDistCum = nResPerDistCum ./ nPerDistCum;
jkl=1;

