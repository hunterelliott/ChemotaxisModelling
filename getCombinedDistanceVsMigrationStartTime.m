function [dC,tC,tCstd,frTrk] = getCombinedDistanceVsMigrationStartTime(t,d,varargin)
%this function returns a single distance vs time curve for the input
%per-cell distance vs. migration initiation time data. It excludes "late
%responders" and bins points by distance
%
% [dC,tC,tCstd,frTrk] = getCombinedDistanceVsMigrationStartTime(t,d,varargin)
%
%   Input:
%       t,d - vectors of migration start times and corresponding distances
%
%   Output:
%
%       dC,tC - vector of combined start times and distances (using
%       average of "first responder" cells)
%
%       dCstd,tCstd - std of first responder times/distances
%
%       frTrk - true/false vector specifying if track was included in
%       first-responders
%
%
%NOTE: All distances are in the same units as d, which isusually microns.
%
%
%Hunter Elliott
%2/2015

%% -------- Input ----------- %%

ip = inputParser;
ip.addParameter('DistanceBinSize',10,@(x)(isscalar(x) && x > 0));%Size to bin distances by
ip.addParameter('MinimumNum',0,@(x)(isscalar(x) && x >= 0));%Minimum number of points in a bin to include in final curve
ip.addParameter('FirstResponderPercentile',10,@(x)(isscalar(x) && x > 0  && x <= 100));%This percentage of the first cells to respond will be considered the "first responders" and included in the curve
ip.addParameter('ShowFigures',true,@islogical);%If true, a figure showing the raw data and extracted curve will be displayed
ip.addParameter('MaximumDistance',[],@(x)(isscalar(x) && x > 0));%Maximum valid distance to include data from - should be the shortest distance between target and image boundary.

ip.parse(varargin{:});

p = ip.Results;


if isempty(p.MaximumDistance)
    warning('combDistVsStart:noMaxDist','No maximum distance was specified! This may produce artifacts in d vs t curve!');
    p.MaximumDistance = max(d);
end

if ~isequal(size(t),size(d))
    error('t and d must be equal sized vectors!')
end


%% -------- Init ----------- %%


distBins = 0:p.DistanceBinSize:p.MaximumDistance;
nBins = numel(distBins)-1;
dC = nan(nBins,1);%Combined distances
tC = nan(nBins,1);%Combined times
tCstd = nan(nBins,1);%STD of combined times
dCstd = nan(nBins,1);%STD of combined distances


%Remove tracks with no start time
noStart = isnan(t);
t(noStart) = [];
d(noStart) = [];

nTrks = numel(t);
frTrk = false(nTrks,1);%Keep track of which tracks were actually used

%% -------- Init ----------- %%

%We break this down by distance to ease definition of "first responders"
for j = 1:nBins
    
    inBin = d >= distBins(j) & d < distBins(j+1);
    
    if nnz(inBin) >= p.MinimumNum
    
        %Get time corresponding to specified percentile for "first
        %responder" cells and find matching cells
        maxT = prctile(t(inBin),p.FirstResponderPercentile);        
        
        fstResp = t <= maxT;%Do it this stupid way so it's easy to keep track of who's included        
        ti = t(fstResp & inBin);
        di = d(fstResp & inBin);
                
        tC(j) = mean(ti);        
        tCstd(j) = std(ti);
        dC(j) = mean(di);%Also record distances...
        dCstd(j) = std(di);
        
        %And log which were counted as first responders
        frTrk(inBin & fstResp) = true;
    
    end    
    
end

if p.ShowFigures
    
    figure        
    plot(d(~frTrk),t(~frTrk),'.k');
    hold on
    plot(d(frTrk),t(frTrk),'.g');    
    errorbar(dC,tC,tCstd)
    legend('Excluded','Included ("First Responders")','Combined Curve +/- STD')
    xlabel('Distance')
    ylabel('Time')
end
    


