function chemotaxisModelFitFile(trackFile,imFile,varargin)
%CHEMOTAXISMODELFITFILE fits diffusion and chemotaxis model to input cellt rack data
%
% chemotaxisModelFitFile(trackFile,imFile)
% chemotaxisModelFitFile(trackFile,imFile,'OptionName1,optionValue1,...)
%
% Input:
%
%   trackFile - name of imaris track position .csv file 
%
%   imFile - name of image file tracks were derived from (for determining
%            target position)
%
%       *See code for this function for details of additional options
%
% Output: saved to disk in same directory as track file unless otherwise
% specified


%% ----------- Input ------------ %%

if nargin < 1
    trackFile = '';
end

[trackPath,trackFile] = optionalFileInput(trackFile,'*.csv','Select an imaris track position file:');


if nargin < 2
    imFile = '';
end

[imPath,imFile] = optionalFileInput(imFile,'*.*','Select the image file:');

ip = inputParser;
ip.addParameter('ReleaseModel','SingleRelease');%Chemoattractant release model type to fit. See fitDiffusionModelToChemotacticPersistence for details.
ip.addParameter('UseMultiStart',true,@islogical);%If true, MultiStart "global optimization" is used
ip.addParameter('NumStart',1e3,@isposint);%Number of initial conditions to use for multistart
ip.addParameter('NumWorkers',6,@isposint);%Number of workers to use for multistart
%ip.addParameter('DistanceBinSize',20,@(x)(isscalar(x) && x > 0));%Size of distance bins to use when creating chemotactic persistence map
ip.addParameter('OutputDirectory',[trackPath trackFile(1:end-4) '_chemoModelling'],@(x)(exist(x,'dir')));%Directory to store analysis results in. If not specified a directory will be created adjacent to the track file
ip.addParameter('TargetChannel',2,@isposint);%Channel number in image file which contains image of target
ip.addParameter('InitialGuess',[],@isstruct);%Optionally input structure with initial guess for modelling parameters
ip.addParameter('BatchMode',false,@islogical);%If true, figure display is supressed for batch processing.
ip.parse(varargin{:});
p = ip.Results;

if ~exist(p.OutputDirectory,'dir')
    mkdir(p.OutputDirectory)
end

if isempty(p.InitialGuess)
    %We have defaults here since our attempted global optimization should
    %make these relatively unimportant (within a couple orders of magnitude
    %or so)
    p.InitialGuess.Cci = 10e-9; %rough estimated values from Daniel
    p.InitialGuess.Di = 1.5e-10;    
    p.InitialGuess.ti = 300;
    p.InitialGuess.Mi = 1e-19;
    p.InitialGuess.cpMaxi = 1;
    p.InitialGuess.cpMini = 0;
    p.InitialGuess.cpSi = 1/p.InitialGuess.Cci;
    p.InitialGuess.Ei = 3/2;
    
end

%% ------ Chemotaxis Analysis ------ %%

%Load the tracks
[trk,frTimes] = imarisReadTrackCSV([trackPath trackFile]);


% --- Find Target Position --- %

[targPos,minValidDist,maxValidDist] = getChemotaxisTargetPosition([imPath imFile],trk,p.TargetChannel,p.OutputDirectory,p.BatchMode);

% --- Run per-track chemotaxis analysis ---- %

trk = trackChemotaxisAnalysis(trk,targPos,'FrameTimes',frTimes);


% --- And get chemotactic persistence map --- %

disp('Making chemotactic persistence map...')
tic
[cpMap,distBins] = migrationParamVsDistance(trk,targPos,frTimes,'ChemoPers','ShowFigures',true,'MaxDistance',maxValidDist,'MinDistance',minValidDist,'BatchMode',p.BatchMode);
mfFigureExport(gcf,[p.OutputDirectory filesep 'chemotactic_persistence_map'])
disp('Finished')
toc


%Save track analysis to output dir
save([p.OutputDirectory filesep 'track_chemotaxis_analysis.mat'],'trk','frTimes','targPos','trackFile','imFile','cpMap','distBins','maxValidDist','minValidDist','p')



%% ----- Model Fitting ---- %%

%Get distance values for model evaluation - we use bin centers
dX = diff(distBins(1:2));
distVals = (minValidDist+dX/2):dX:(maxValidDist-dX/2);


fitInfo = fitDiffusionModelToChemotacticPersistence(frTimes,distVals,cpMap,p.InitialGuess.Di,p.InitialGuess.Mi,p.InitialGuess.Cci,... %TEMP - Ugh. This call is hideous - convert to input structure passing?
    p.InitialGuess.ti,p.InitialGuess.cpMaxi,p.InitialGuess.cpMini,p.InitialGuess.cpSi,'Ei',p.InitialGuess.Ei,'ReleaseModel',...
    p.ReleaseModel,'UseMultiStart',p.UseMultiStart,'NumStart',p.NumStart,'NumWorkers',p.NumWorkers); %#ok<NASGU>


outFile = [p.OutputDirectory filesep 'diffusion_model_fit_' p.ReleaseModel];
if p.UseMultiStart
    outFile = [outFile '_MultiStart' num2str(p.NumStart)];
end

save([outFile '.mat'],'fitInfo','cpMap','distVals','distBins','frTimes','p','trackFile','imFile','maxValidDist','minValidDist','p')



