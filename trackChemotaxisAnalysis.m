function trk = trackChemotaxisAnalysis(trk,tPos,varargin)

% Input:
% trk- structure with 2D track info, in format produced by e.g. imarisReadTrackCSV
% tPos - 2 element vector with chemotaxis target position
%
%
% Output:
%
%   trk - track structure with additional per-track chemotaxis parameters calculated

%Hunter Elliott
%?/2015


ip = inputParser;
ip.addParameter('FrameTimes',[]); %Times in seconds corresponding to each frame in the movie the tracks were derived from
ip.addParameter('ShowFigures',false,@islogical); %Show figures with track parameters
ip.parse(varargin{:});
p = ip.Results;

nTrk = numel(trk);

if p.ShowFigures
    allFig = figure;
    hold on
end

if isempty(p.FrameTimes)
    %Just use all times present in tracks
    p.FrameTimes = unique(vertcat(trk.Time));    
    warning('Frame times not input, assuming all timepoints are represented in tracks!')
end


%% ----- Parameters ----= %%

dTol = 1;
minPts = 10;%Minimum points per track. We're extracting a shit ton of parameters from each so need a few points....

dDistFit = cell(nTrk,1);
% fitOpt = fitoptions('smoothingspline');
% fitOpt.SmoothingParam = 1e-5;
distTol = 2.5;%Roughly a cell radius
dDistTol = -.075;%Arbitrary for now - set by an example non-directional cell.
minRunLen = 60;%Using moving average size from cited paper.
tRes = 1;%Sub-frame resolution to use
minPers = 2/3;%Minimum chemotactic persistence to be considered directional


%% -------- Per-track chemotaxis analysis ------ %%

aDxDtAll = nan(nTrk,max(vertcat(trk.Frame)));
aDxDtAllMA = nan(nTrk,max(vertcat(trk.Frame)));
distTAll =  nan(nTrk,max(vertcat(trk.Frame)));


dStart = nan(nTrk,1);
nPtsPer = nan(nTrk,1);

tTot = tic;
disp('Running per-track chemotaxis analysis...')

for j = 1:nTrk
    
    nPtsPer(j) = numel(trk(j).x);
    
    if nPtsPer(j) > minPts
    
        % ---- Rate of advance towards target calc ---- %    
        X = [trk(j).x,trk(j).y];
        t = trk(j).Time;
        dx = diff(X,1,1);%Displacement
        dxN = bsxfun(@times,dx,1 ./ sqrt(sum(dx .^2,2)));%Normalized direction vector
        dt = bsxfun(@minus,tPos,X);%Vector towards target
        distT = sqrt(sum(dt .^2,2));%Distance to target
        dtN = bsxfun(@times,dt,1 ./ distT);%Normalized target direction vector

        %Get rate of change of distance to target
        dDistT = diff(distT) ./ diff(t);                

        %Fit to the positions so we can more easily set a tolerance
        currTol = numel(dDistT)*distTol^2/2;%Get tolerance for current fit TEMP - this is rough currently, designed to get 95% of errors to be < cell radius
        distTFit = spaps(t,distT,currTol);
        dDistTFit = fnder(distTFit,1);        

        trk(j).dDistToTargetSmooth = fnval(dDistTFit,t);
        trk(j).DistToTarget = distT;
        trk(j).DistToTargetSmooth = fnval(distTFit,t);

        tSubRes = min(t):tRes:max(t);
        dDistSubRes = fnval(dDistTFit,tSubRes);

        isAdv = dDistSubRes<dDistTol;

        %We always preserve events that are truncated by dataset end
        if isAdv(end)        
            isAdvL = bwlabel(isAdv);
            endAdv = isAdvL == isAdvL(end);
        else
            endAdv = false(size(isAdv));
        end    
        isAdv = bwareaopen(isAdv,round(minRunLen/tRes)) | endAdv;



        % ----- Chemotactic Persistence ----- %

        cumDisp = vertcat(0,cumsum(sqrt(sum(dx .^2,2))));%Cumulative displacement
        cumDispFit = spaps(t,cumDisp,currTol);
        cumV = fnder(cumDispFit,1);%Smoothed velocity        
        cumVSubRes = fnval(cumV,tSubRes);

        chemPers = -dDistSubRes ./ cumVSubRes;
        
        %Store at original time res for later display/analysis
        trk(j).ChemoPers = -trk(j).dDistToTargetSmooth ./ fnval(cumV,t);        
        trk(j).v = sqrt(sum(dx .^2,2)) ./ diff(t);
        trk(j).smoothV = fnval(cumV,t);

        isPers = chemPers > minPers;

        % ----- Relative Angle ----- %

        %Angle between displacement and target separation vector
        trk(j).ChemoInd = dot(dtN(1:end-1,:),dxN,2);%Cos of angle between velocity and direction to target (chemotactic index)
        aDxDt = acos(trk(j).ChemoInd);%Angle between displacement and target separation vector
        aDxDtMA = smooth(aDxDt,5);

        currFrames = trk(j).Frame;
        aDxDtAll(j,currFrames(1:end-1)) = aDxDt;
        aDxDtAllMA(j,currFrames(1:end-1)) = aDxDtMA;
        distTAll(j,currFrames) = distT;

        % ----- Start point determination ----- %m

        if nnz(isAdv & isPers) > 0
            iStartSubRes = find(isAdv & isPers,1,'first');
            trk(j).ChemoStartT = tSubRes(iStartSubRes);
            dStart(j) = fnval(distTFit,trk(j).ChemoStartT);
            trk(j).ChemoStartIPer = min(max(round(iStartSubRes / numel(tSubRes) * numel(t)),1),numel(t));%The per-track index of start.
            [~,trk(j).ChemoStartI] = min(abs(p.FrameTimes-trk(j).ChemoStartT));        
            %trk(j).ChemoStartT = t(trk(j).ChemoStartI);
            %dStart(j) = distT(trk(j).ChemoStartI);
        end
        
         

        if p.ShowFigures
            figure
            hold on

            plot(trk(j).x,trk(j).y,'k')
            scatter(trk(j).x,trk(j).y,5,trk(j).Time)
            axis equal

            plot(tPos(1),tPos(2),'rx')

            quiver(X(:,1),X(:,2),dtN(:,1),dtN(:,2),0)       

            figure
            hold on

            plot(trk(j).x,trk(j).y,'k')
            scatter(trk(j).x(1:end-1),trk(j).y(1:end-1),50,trk(j).ChemoInd)
            axis equal

            plot(tPos(1),tPos(2),'rx')

            quiver(X(:,1),X(:,2),dtN(:,1),dtN(:,2),0)       
            
            
            figure
            plot(aDxDt,'.-')
            ylabel('Angle to target')
            hold on
            plot(aDxDtMA,'.-g')

            figure
            subplot(4,1,1)        
            plot(t,distT,'.-r')
            ylabel('Distance to target')
            hold on
            fnplt(distTFit)
            subplot(4,1,2)
            plot(t(1:end-1),dDistT,'.-r')
            ylabel('delta distance to target')
            hold on
            plot(xlim,[0 0 ],'--k')
            plot(xlim,[1 1]*dDistTol,'--r')
            fnplt(dDistTFit)
            %plot(dDistFit{j},'.-r')
            %plot(aDxDtMA,'.-g')
            subplot(4,1,3)
            plot(t,cumDisp,'r.-')
            hold on
            fnplt(cumDispFit)
            ylabel('Cumulative path length')

            subplot(4,1,4)
            plot(tSubRes,chemPers,'b.-')
            hold on
            plot(t(1:end-1),cumDisp(2:end) ./ dDistT,'.-r')
            ylabel('Ratio, delta path length to delta d-to-target')
            ylim([-1 1])
            plot(xlim,[0 0],'--k')
            plot(xlim,[1 1] * minPers,'--r')
            xlabel('Time, seconds')

        end
    end        
    
end

disp(['Finished chemotaxis analysis! Took ' num2str(toc(tTot)) ' seconds'])

