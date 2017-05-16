function fitInfo = fitDiffusionModelToChemotacticPersistence(t,d,cp,Di,Mi,Cci,tii,cpMaxi,cpMini,cpSi,varargin)

ip = inputParser;
ip.addParameter('dT',1,@(x)(isscalar(x) && x > 0));%Time resolution to use for non-discrete release models.
ip.addParameter('ReleaseModel','GeometricScaling',@(x)(ismember(x,{'GeometricScaling','SingleRelease','ContinuousRelease','PulseRelease'})));
ip.addParameter('Ei',3/2,@(x)(isscalar(x) && x >= 0));%Initial value of exponent in geometric scaling model
ip.addParameter('te',tii+60,@(x)(isscalar(x) && x >= 0));%Initial value of pulse end time in PulseRelease model
ip.addParameter('UseMultiStart',true,@islogical);%If true, use multi-start to attempt to find global optimal fit
ip.addParameter('NumStart',50,@(x)(isposint(x) && x > 1));%Number of start points to use if MultiStart enabled
ip.addParameter('NumWorkers',6,@(x)(isposint(x)));%Number of parallel workers to use if multi-start enabled
ip.addParameter('Solver','lsqnonlin',@(x)(ismember(x,{'lsqnonlin'})));%Choose solver to use for fit. Not used currently, just future-proofing - Current objective function only supports lsqnonlin
ip.parse(varargin{:});
paramIn = ip.Results;


%% ---- Initialization ---- %%

%Initial value vector
pI = [Di, Mi, Cci, tii, cpMaxi, cpMini,cpSi];

%Define bounds
lb = [Di / 1e2, Mi/1e2, Cci / 1e2 , 0 , 1e-2, 0, 1e-3/Cci];
ub = [Di * 1e2, Mi*1e2, Cci * 1e2 , max(t)-diff(t(1:2)), 1, .9, 1e3/Cci];
% lb = [];
% ub = [];

%Scale parameters to similar orders of magnitude to avoid numerical problems in optimization
pScaling = [1e10 1e18 1e8 1e-2 1 1 1e-8];

switch paramIn.ReleaseModel
    
    case 'GeometricScaling'
        
        %In this case we fit the exponent as well.        
        extraPar = {'E'};                
        pI = [pI paramIn.Ei];
        lb = [lb 0];
        ub = [ub paramIn.Ei*10];
        pScaling = [pScaling 1];

    case 'PulseRelease'
        %This model needs to fit pulse end time as well        
        extraPar = {'te'};
        pI = [pI paramIn.te];
        lb = [lb min(t)+paramIn.dT];
        ub = [ub max(t)];
        pScaling = [pScaling 1e-2];

    otherwise
        extraPar = {};
        
end

 
opts = optimoptions('lsqnonlin');
%opts.Algorithm = 'levenberg-marquardt';

% opts.DiffMaxChange = 1e-9;
%opts.TypicalX = pI;
opts.Diagnostics  = 'on';
if paramIn.UseMultiStart    
   opts.Display = 'off';
else
    opts.Display = 'iter-detailed';
end
%opts.TolFun = 1e-7;
opts.FinDiffType = 'central';
%opts.FinDiffRelStep = 1e-15;
%opts.TolX = 1e-200;
opts.MaxFunEvals = 5e3;
opts.MaxIter = 5e2;

%% ---- Fitting ---- %%

fitFun = @(p)(fitFunction(p,t,d,cp,pScaling,paramIn.ReleaseModel,paramIn.dT));

disp(['Fitting ' paramIn.ReleaseModel ' diffusion model...'])

optProb = createOptimProblem(paramIn.Solver,'objective',fitFun,'lb',lb .* pScaling,'ub',ub .* pScaling,'x0',pI .* pScaling,'options',opts);

tic
if paramIn.UseMultiStart
    %Use MultiStart to attempt to find global optimal solution.
          
    ms = MultiStart('UseParallel',paramIn.NumWorkers > 1,'Display','iter');
    
    if paramIn.NumWorkers > 1
        %Setup pool if necessary
        poolObj = gcp('nocreate');
        if isempty(poolObj)
            %If it already exists we just use that size.
            poolObj = parpool(paramIn.NumWorkers)
        end
        
    end
    
    
    [pFit,~,~,outputMS,sltnsMS] = run(ms,optProb,paramIn.NumStart);
    
    %Store the info about the MultiStart run
    fitInfo.AllSolutions = sltnsMS;
    fitInfo.OutputMultiStart = outputMS;
    
    %Now use the best solution to run lsqnonlin for approximation of jacobian at best solution
    optProb.x0 = pFit;
    
    
end
    

[pFit,resN,resid,exFlag,output,lambda,J] = lsqnonlin(optProb);

  
    
fitInfo.Res = resid;        
fitInfo.Lambda = lambda;
fitInfo.J = full(J);
fitInfo.ResNorm = resN;
fitInfo.exFlag = exFlag;
fitInfo.Output = output;

disp('Done.')
toc

%% ---- Parameter covariance and CI estimates ---- %%


% Calculate covariance matrix
[~,R] = qr(fitInfo.J,0);
Rinv = R\eye(size(R));
pInvJtJ = Rinv*Rinv';
fitInfo.Cov = var(fitInfo.Res(:)) * pInvJtJ;
fitInfo.Cov = fitInfo.Cov ./ (pScaling' * pScaling);%Undo scaling

%undo this scaling last to avoid numerical issues
fitInfo.J = bsxfun(@rdivide,fitInfo.J,pScaling);

try
    [std,corr] = cov2corr(fitInfo.Cov);
    fitInfo.Std = std;
    fitInfo.Corr = corr;
catch %TEMP - avoid singular matrix return

end

%Parameter confidence intervals
pCi = nlparci(pFit ./ pScaling,resN,'covar',fitInfo.Cov);
    

%Undo scaling, store in output for posterity
pFit = pFit ./ pScaling;
fitInfo.ParScaling = pScaling;%Useful in interpreting multistart output.


%Put it all in a structure for output
fitInfo.ParOrder = [{'D','M','Cc','ti','cpMax','cpMin','cpS'} extraPar];

for j = 1:numel(pFit)
    fitInfo.(fitInfo.ParOrder{j}) = pFit(j);    
    fitInfo.([fitInfo.ParOrder{j} '_CI']) = pCi(j,:);    
end

%Also return fitted concentration, persistence and associated distances and
%times
[~,fitInfo.ConcFit,fitInfo.ChemPersFit] = fitFun(pFit .* pScaling);

%We need to calculate residuals ourselves for MultiStart
if paramIn.UseMultiStart
    fitInfo.Res = cp - fitInfo.ChemPersFit;
end



end

function [res,c,cpE] = fitFunction(p,t,d,cp,pScaling,releaseModel,dT)        
        
    p = p ./ pScaling;
    
    switch releaseModel
        
        case 'GeometricScaling'
            %The parameter only specifies the beginning of release and rate
            %of mass release in this case.
            ti = p(4):dT:max(t);
            mi = p(2) .* (ti / ti(1)) .^ p(end);
            
        case 'SingleRelease'
            
            ti = p(4);
            mi = p(2);
            
        case 'ContinuousRelease'
            
            ti = p(4):dT:max(t);
            mi = ones(1,numel(ti)) * p(2);
            
        case 'PulseRelease'
            
            ti = p(4):dT:p(end);
            mi = ones(1,numel(ti)) * p(2);
            
        otherwise
            error('unrecognized release model name')
    end
    %Get simulated concentrations for these parameters
    c = diffusionModelRadial3DPlanarBoundPointSource(t,d * 1e-6,p(1),ti,mi);
    %And expected corresponding chemotactic persistence using sigmoidal
    %response model
    cpE = p(6) + (p(5) - p(6)) ./ (1+exp(-p(7) .* (c - p(3))));
    
    %And residuals
    res = cp - cpE;
            

end


