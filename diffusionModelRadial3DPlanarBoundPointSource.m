function c = diffusionModelRadial3DPlanarBoundPointSource(t,r,D,ti,mi)
%DIFFUSIONMODELRADIAL3DPLANARBOUNDPOINTSOURCE models radial concentration distribution after point release(s) of diffusive mass
%
% c = diffusionModelRadial3DPlanarBoundPointSource(t,r,D,ti,mi)
%
%   This model assumes diffusion within a 3D domain with a planar boundary
%   at z=0, at time t after one or more point-source releases of diffusive
%   material.
%
% Inputs:
%
%   t - scalar  or vector of times to calculate c at
%   r - scalar or vector of radii to calculate c at 
%   D - diffusion constant
%   ti - scalar or vector of time(s) of source release
%   mi - scalar or vector of mass(es) released at times ti.If scalar and ti
%   is vector an equal mass at each ti is assumed.
%
% Output:
%
%   c - RxT matrix with the concentration at all R radii and all T times
%

%Hunter Elliott
%2/2015


nTi = numel(ti);%Number of releases
nR = numel(r);%Number of radii
nT = numel(t);%Number of times

if isscalar(mi) && nTi > 1
    %Assume equal mass released each time
    mi = repmat(mi,[nTi, 1]);
end

%Define model for a single release at tr
cFun = @(t,tr,r,mi)( (2 .* mi ) ./ (4*pi*D .* (t - tr)) .^ (3/2) .* exp(-r .^2 ./ (4*D .* (t-tr))));

%Make matrices for vectorization
rMat = repmat(r(:)',[nTi, 1, nT]);
tMat = repmat(reshape(t,[1 1 nT]),[nTi nR 1]);
tiMat = repmat(ti(:),[1 nR, nT]);
miMat = repmat(mi(:),[1 nR, nT]);

%Get individual release concentrations
cPer = cFun(tMat,tiMat,rMat,miMat);
%Handle times prior to and at release as a special case. Because t=tr is a singularity we also exclude it here
cPer(tiMat >= tMat) = 0;
%And final solution using principle of superposition
c = squeeze(sum(cPer,1));

