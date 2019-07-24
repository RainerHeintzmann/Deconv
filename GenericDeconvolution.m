% [res,resIllu,resPSF,evolObj,evolIllu,evolOTF,allObj,allIllu,allOTF]=GenericDeconvolution(image,psf,NumIter,Method,Update,Regularisation,betas,borderSizes,Variances,useCuda,keepExtendedStack) : Deconvolution routine with various functionalities.
% image: Image to Deconvolve
% psf: Point spread function to deconvolve with
% NumIter(1): Iterations to perform. Optionally the subiterations for initial Object iterations (2) and successive object (3), illumination (4) and psf (5) iterations  can be given.
% Method: One of 'LeastSqr', 'WeightedLeastSqr', 'Poisson' or 'GaussianWithReadnoise', 'Empty' defining the norm of the data-simulation agreement to optimize. 'Empty' always returns zero (useful for gradient tests of the regularizers)
% 
% Update: Default: [] uses 'lbfgs', This describes the minimization algorithm (update scheme) to use. You can use all the ones that minFunc allows, but in addition also 'RL' for
% RichardsonLucy multiplicative (EM) Algorithm, or 'RLL' for a Wolfe line search along the Richardson Lucy update direction
% 'RLR' uses the recursive way of accelerating the RL algorithm based on Biggs&Andrews, Applied Optics 36,1-9.
% Regularization : A cell array containing the regularisers to use and their respective arguments. One list can be specified for each of object, illumination and psf e.g. {{},{},{}} 
% E.g. {'TV',[0.02 1.0]} will use total variation with a lambda of 0.02 and an epsilon of 1.0
%       They can be combined (e.g. {'TV',[0.02 1.0];'NegSqr',0.05}. 
%       In the case of blind deconvolution two such respective arrays have to be provided. E.g. {{'TV',[0.02 1.0];'NegSqr',0.05},{'NegSqr'}}.
% In this Regularization list also the following flags are possible:
%       'TV',[lambda, epsilon]: Total variation modified according to Ferreol Soulez et al. Blind deconvolution of 3D data in wide field fluorescence microscopy
%                               choose epsilon=0 for standard TV. penalty = sqrt(|grad(f)|^2+epsR)
%       'GS',[lambda]: Gradient squared regularisation. penalty= |Grad(f)|^2
%       'GR',lambda: Good's roughness regularisation. penalty= |Grad(f)|^2 / (|f|+epsR)
%       'CO',lambda: Conchello's |f|^2 regularisation (Gaussian prior)
%       'RealSpaceMask' used for light sheet microscopy in the PSF regularisation field. If only a number is supplied, this is interpreted as the width of the Gaussian
%       'ER',[lambda,gamma]: ER-regularization according to Arigovindan et al., PNAS, 2013. penalty=ln(epsR+ (|f|^2+ SumHessianSqr)),  SumHessianSqr= (d^2/dxdx f)^2 + (d^2/dydy f)^2+ (2d^2/dxdy f)^2
%       'Bg', value: Fixed background value (to get the Poisson Statistics right)
%       'ForcePos',[] : Will make the object the square of an auxiliary function, which is then estimated
%       'ForcePhase',[] : Will make the object exp(1i*f) of an auxiliary function f, which is then estimated
%       'Complex',[] : Will allow the object (if an auxiliary function is used, the auxiliary function) to be an aplitude object. Note that even for complex valued data the object can remain to be real.
%       'Resample',[sx sy sz] : Specifies the resampling to be used for reconstructed object (only for object at the moment)
%       'NormMeasSum',[] : Normalizes the result of the FwdModel to the sum of all measured values befor calculating the error and residuum
%       'NormMeasSumSqr',[] : Normalizes the result of the FwdModel to the sqrt of the sum of squares of all measured values befor calculating the error and residuum
%       'NormFac',[NF] : Allows to set the NormFactors (rather than automatic determination. Negative values mean automatic determination for the respective estimation steps.
%       'Reuse' : The algorithm will use a previously caculated result (stored in the global variable aRecon) and continue from there
%       'StartImg', anImage : Uses this image as a starting object for the iterations
%       'Illumination', anIllumination : Uses the illuminations as given 
%       'IlluMask', aMask : Allows the user to provide an illumination mask
%       'MaxTestDim',  for the Gradien testing (use negative iterations number, which amounts to the epsilon for the numeric steps!). How many directions to test in real before going to imaginary.
%
% betas: Vector of beta values to vary the weight of different directions (e.g. anisotropic sampling). This should correspond to pixelsizes, for example normalized to the x-pixel size
% or optical resolution
% borderSizes: If given, the reconstruction is done over an extended area, but the non-existing data there is assumed to be always correct (useful especially for widefield data)
%              Note that you can also give a boolean map here, which notes which of the data is valid.
% Variances: this is useful to give a map of expected variance (when using the 'LeastSqr' method)
% useCuda: This can accelerate the reconstruction significantly (the CudaMat toolbox by Rainer Heintzmann has to be installed)
%
% Requires: DipImage and Mark Schmid's minFunc (http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html).
% It can work with CudaMat to accelerate the deconvolution 
% In the file polyinterp.m line 103 needs to be changed to:
% for qq = 1:length(cp); xCP = cp(qq);
% If everything is prepared run:
% speedtestDeconv(0)
% e.g. result could be 29.4 sec and
% speedtestDeconv(1)
% e.g. result could be 3.11 sec
%
% Example calls:
%
% %2D Example for Poisson noise, Goods roughness, positivity constraint and 10 pixels boardser:
% obj=readim; MaxPhotons=1000;cutSize=[200,200];obj=obj/max(obj)*MaxPhotons;
% psf=abssqr(ift(rr(size(obj))<30)); otf=ft(psf);otf=otf/max(abs(otf));
% img=extract(noise(ift(ft(obj) .* otf),'poisson'),cutSize); useCuda=0;
% myDeconvGP=GenericDeconvolution(img,psf,35,'Poisson',[],{'GR',[0.005 0.1];'ForcePos',[];},[1,1],[10 10],[],useCuda); gtp=cat(1,img{1},myDeconvGP,extract(obj,cutSize))
% myDeconvTV=GenericDeconvolution(img,psf,35,'Poisson',[],{'TV',[0.0004 500];'ForcePos',[];},[1,1],[10 10],[],useCuda); gtp=cat(1,img{1},myDeconvTV,extract(obj,cutSize))
% Example use for blind PSF deconvolution:
%    [obj,dum,psf]=GenericDeconvolution(img,psf,[2 0 40 0 70],'LeastSqr','Blbfgs',{{'StartImg',obj;'ForcePos',[]},{},{'NegSqr',1e9;'CO',1e7}},[1,1,1],[0 0 0],[],useCuda); 
% %For complex data and complex PSF
% myDeconv=GenericDeconvolution(img,h,15,'LeastSqr',[],{'Complex',[]},[1 1 1],[0 0 0],[],useCuda); st=cat(1,img{1},myDeconv)
%
% %An example with 2x oversampling in the reconstruction:
% myDeconvGP=GenericDeconvolution(img,h,85,'Poisson',[],{'GR',[0.01 0.1];'ForcePos',[];'Resample',2},[1,1,1],[0 0 0],[],useCuda); gtp=cat(1,img{1},myDeconvGP)
%
% Example for a gradient test:
% obj=readim; MaxPhotons=100;cutSize=[20,20];obj=obj/max(obj)*MaxPhotons;
% psf=abssqr(ift(rr(size(obj))<5)); otf=ft(psf);otf=otf/max(abs(otf));
% img=extract(noise(abs(ift(ft(obj) .* otf)),'poisson'),cutSize); useCuda=0;
% myDeconv=GenericDeconvolution(img,h,-0.1,'Empty',[],{'MaxTestDim',10;'NormFac',1;'GR',0.1;'StartImg',img},[1,1,1],[0 0 0],[],useCuda);%
%
% %Contributing Authors: Rainer Heintzmann, Aurelie Jost, Polina Feldmann


%***************************************************************************
%   Copyright (C) 2008-2009 by Rainer Heintzmann                          *
%   heintzmann@gmail.com                                                  *
%                                                                         *
%   This program is free software; you can redistribute it and/or modify  *
%   it under the terms of the GNU General Public License as published by  *
%   the Free Software Foundation; Version 2 of the License.               *
%                                                                         *
%   This program is distributed in the hope that it will be useful,       *
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
%   GNU General Public License for more details.                          *
%                                                                         *
%   You should have received a copy of the GNU General Public License     *
%   along with this program; if not, write to the                         *
%   Free Software Foundation, Inc.,                                       *
%   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
%**************************************************************************
%

function [res,resIllu,resPSF,evolObj,evolIllu,evolOTF,allRObj,allRIllu,allROTF]=GenericDeconvolution(image,psf,NumIter,Method,Update,Regularisation,betas,borderSizes,Variances,useCuda,keepExtendedStack)
global myim; % Cell array of the images
global myFTim; % used only for the Ptychography trick to avoid the back transforms
global otfrep; % Cell array of the otfs
global DeconvMethod;
global aResampling;   % before convolution
global subSampling;   % after convolution 
global RegularisationParameters;
global DeconvMask;
global DeconvVariances;
global BetaVals;
global myillu; 
global myillu_mask; 
global OTFmask;
global ToEstimate;   % This flag controls what to estimate in this iteration step. 0: sample density, 1: illumination intensity, 2: both, 3: psf
global aRecon;   % During the deconvolution the object is stored as a dipImage in this variable
global NormFac;  % normalisation factor
global myillu_sumcond;
% global TheObject;
global cuda_enabled;  % also to see if cuda is installed
global useCudaGlobal;
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF
global ConvertModelToVec; % Contains the routine to convert the model (e.g. the gradient of the object) into the vector to iterate
global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
global RefObject;  % If supplied by the user, it is used to check the progress during iterations
global RefObject_SSQ;  % For checking the progress to the ground truth during iterations.  RH 2017
global RefObject_SAbs; % For checking the progress to the ground truth during iterations.  RH 2017
global RefObject_SSQ2;  % For checking the progress to the ground truth during iterations.  RH 2017
global RefObject_SAbs2; % For checking the progress to the ground truth during iterations.  RH 2017

%global ComplexObj;    % Can the object to estimate be complex valued?
%global IntensityData;  % Is the supplied data the abs square of the object to estimate?
%global ForcePos;
%ForcePos=0;         % Will only be set to one in the ParseRegularisation routine, if chosen.
global ComplexPSF; % will be set, if PSF is complex valued
global PupilInterpolators;  % used for estimating ATFs in blind deconvolution
global measSumsSqr;
global measSums;
global allObj;
global ReadVariance;
ReadVariance=1.0;  % Just to give it a default value, if the user does not use the parameter "ReadVariance"

myFTim=[];
DeconvMethod=Method;
aResampling=1; 
subSampling=1; 
%if ~isempty(PupilInterpolators) 
%end
msevalue=0;
evolObj=[];
evolIllu=[];
evolOTF=[];
if nargout > 6
    allObj = {};
else
    allObj = [];
end

if nargin < 11
    keepExtendedStack=0;
end
if nargin < 10
    useCuda=0;
end
useCudaGlobal=useCuda;

if nargin < 9
    Variances=[];
end
if nargin < 8
    borderSizes=[0 0 0];  % These can help to avoid artefacts
end
if nargin < 7
    betas=[1 1 1];
end
if nargin < 6
    Regularisation={};  % No regularisation performed
end
if nargin < 5 || isempty(Update)
    Update='lbfgs';   % Limited memory BFGS
end
if nargin < 4
    Method = 'Poisson';
end
if nargin < 3
    NumIter=50;
end

if size(NumIter) < 2
    NumIter(2) = 5;  % Default value for SIM iterations. For Speckles use 0
end

if size(NumIter) < 3
    NumIter(3) = 25;  % Default value for SIM iterations. For Speckles use 5
end

if size(NumIter) < 4
    NumIter(4) = 5;  % Default value for SIM iterations. For Speckles use 15
end

if Update(1)=='B'
    myillu=[];  % force the illumination to be empty. Also for the correct normalizations of start vectors
end

myim=image;  % Just temporarily to be used for size determination
% if isempty(myillu)
%     AssignFunctions(RegularisationParameters,0); % To estimate is Object (non-blind)
% else
%     AssignFunctions(RegularisationParameters,1); % To estimate is Object but fixed illumination has to be considered
% end

%%
% somehow the line below gets overwritten, so we have to set the boudary
% in the inner loop "GenericErrorAndDeriv"
dip_setboundary('periodic')  % Establishes Periodic Boudary conditions.

if isa(borderSizes,'dip_image')
    fprintf('Mask provided as bordersize. Using mask instead of generating borders');
    DeconvMask=borderSizes;
    borderSizes=0 .* size(DeconvMask);
else
    borderSizes=borderSizes*2;
end

if useCuda
    % enableCuda();
    initCuda();   % set all cuda use-variables to 1
    set_ones_cuda(useCuda); set_zeros_cuda(useCuda);   % make sure that the zeros and ones function from cuda is used.
else
    if ~isempty(cuda_enabled)
        disableCuda();  % just to make sure that newim and so on do not generate cuda objects
    end
end

    s=size(image{1});
    if norm(borderSizes) <= 0
        borderSizes=s*0;
    end
    if length(borderSizes) > length(s)
        s(end+1:length(borderSizes))=1;
    end

    extraBorder(1)=0;
    if numel(borderSizes) > 1
        extraBorder(2) = mod(s(2)+borderSizes(2),2);  % Only expand for uneven sizes along Y
    end
    if numel(borderSizes)> 2 % && borderSizes(3)==0
        extraBorder(3)=0;
    end
    if norm(extraBorder) > 0
        fprintf('Warning: Bodersize needed to be extended, to yield and even sized (along Y) total array, required for rft.\n')
        borderSizes=borderSizes+extraBorder;
    end
    clear extraBorder;
    
%myillu={};
NumViews=1;
% lambdaPenalty=lambda;
%RegularisationMethod=Prior;
%NegPenalty=NegPrior;
DeconvMask=[];
DeconvVariances=Variances;
BetaVals=betas;
ToEstimate=0;
NormFac=1;

myim={};  % in case some previous function stored more values in here, we need to clear it.

% unix('touch cudaArith.cu; touch cuda_cuda.c; make'); clear classes  % compile command, just to copy and paste it
if iscell(image) || isa(image,'dip_image_array')
    myim=image;
else
    for v=1:size(image,4)
        if length(size(image))>3
            myim{v}=image(:,:,:,v-1);
        else
            myim{v}=image;
        end
    end
end
clear image;
myim=ConditionalCudaConvert(myim,useCuda); %Seems like there was a case mismatch. Aurelie 03.03.2014

if ~(iscell(psf) || isa(psf,'dip_image_array'))
    for v=1:size(psf,4)
        if length(size(psf))>3
            n{v}=psf(:,:,:,v-1);
        else
            n{v}=psf;
        end
    end
    psf=n;
    clear n;
end
psf=ConditionalCudaConvert(psf,useCuda); %Seems like there was a case mismatch. Aurelie 03.03.2014
OrigSize=size(myim{1});
if length(borderSizes) > length(OrigSize)
   OrigSize(end+1:length(borderSizes))=1;
end
%% Parse the regularisation parameters.  Needs to be before aResampling is used below
[RegObj,RegIllu,RegOTF]=ParseRegularisation(Regularisation);  % Also sets global variable such as ForcePos and ComplexObj
RegularisationParameters=RegObj;

%% Extend the border
NewSize=floor(OrigSize+borderSizes); % Will also be needed for object starting vector
NewDataSize=NewSize;  % Will also be needed for object starting vector
if any(aResampling~=1)
   NewIlluSize = floor(NewSize .* aResampling);
   NewObjSize = floor(NewSize .* aResampling);
else
   NewIlluSize=NewSize;
   NewObjSize=NewSize;
end

if any(subSampling~=1)
   NewIlluSize = NewIlluSize .* subSampling;  % Only integer factors are allowed
   NewObjSize = NewObjSize .* subSampling;    
   NewSize = NewSize .* subSampling;  % for PSF to have the correct size
end

if norm(borderSizes) > 0 || any(subSampling~=1)
        if (1)  % Here the "one-slice speedup" trick can be disabled or enabled
            if length(NewDataSize) > 2 && OrigSize(3) == 1
                fprintf('Border regions: Thick slice speedup possible! Will keep Z-datasize of only one slice.\n');
                NewDataSize(3)=1;
            end
        end
        if norm(OrigSize-NewDataSize) > 0 || (norm(OrigSize-NewSize))
            fprintf('Border region requested, expanding data size from [');
            fprintf('%g,',OrigSize);
            fprintf('] to [');
            fprintf('%g,',NewDataSize);
            fprintf('] \n');
            DeconvMask=extract(myim{1}*0>-1,NewDataSize);  % Place zero outside and one at the old data region
            for v=1:length(myim)  % insert empty information
                myim{v}=extract(myim{v},NewDataSize);  % Overwrite the old data
                if ~isempty(psf{1}) && (v<=numel(psf)) && norm(size(psf{v})-NewSize)~=0  % Not NewObjSize since the Resampling does not affect the OTF
                    psf{v}=squeeze(extract(psf{v},NewSize,[],'cyclic'));  % Not NewObjSize since the Resampling does not affect the OTF
                    fprintf('Warning: View %d, PSF has incorrect size. Adapting the size of the PSF\n',v);
                end
            end
        else
            DeconvMask=[];  % No mask is needed
        end
        if exist('myillu') && ~isempty(myillu)   % illu has to use Object size = NewSize, not NewDataSize
            for v=1:length(myillu)
                sillu=size(myillu{1+mod(v-1,numel(myillu))});
                sillu(end+1:length(NewIlluSize))=1;
                if (v <= numel(myillu)) && norm(sillu-NewIlluSize) ~= 0
                    fprintf('Warning: View %d, The given illumination distribution does not correspond in size to the extended border region. Cyclically expanding to fix this.\n',v)
                    myillu{v}=extract(myillu{v},NewIlluSize,[],'cyclic');
                end
            end
        end
        if ~isempty(Variances)
            DeconvVariances=extract(Variances,NewDataSize,floor(NewDataSize/2),1);  % Pad with Variance of one (gets ignored)
        end
else
    for v=1:numel(psf)   % PSF has to use NewSize, not NewDataSize
        if ~isempty(psf{v})
            spsf=size(psf{1+mod(v-1,numel(psf))});
            OrigSize(end+1:length(spsf))=1;
            if norm(spsf-NewDataSize)~=0
                psf{1+mod(v-1,numel(psf))}=squeeze(extract(psf{v},NewDataSize));   % In case the PSF has a different size
                fprintf('Warning: View %d, PSF has incorrect size. Adapting the size of the PSF\n',v);
            end
        end
    end
end

if (~isempty(PupilInterpolators) && ndims(myim{1})>2 && size(myim{1},3)>1)
    Pupil3DPrepare(size(myim{1}))  % Needs the imagesize to already account for resampling and borders
    
    if ~isempty(PupilInterpolators) && ~isempty(PupilInterpolators.Newlambda)&& ~isreal(psf{1}) %Aurelie 22.10.2014: PupilInterpolators used only in the blind ASF deconvolution case
        PupilInterpolators.indexList2D=ConditionalCudaConvert(PupilInterpolators.indexList2D,useCuda);
        PupilInterpolators.fullIndex3D=ConditionalCudaConvert(PupilInterpolators.fullIndex3D,useCuda);
        PupilInterpolators.factorList=ConditionalCudaConvert(PupilInterpolators.factorList,useCuda);
        PupilInterpolators.Mask=ConditionalCudaConvert(PupilInterpolators.Mask,useCuda,1);
    end
end
%%
if isempty(psf{1}) && ~RegObj(21,1)  % Only FTData does not need PSF
    if RegOTF(14,1)
        ImageParam=struct('Sampling',PupilInterpolators.pixelsize,'Size',size(myim{1}));
        PSFParam=struct('NA',PupilInterpolators.NA,'n',PupilInterpolators.RI,'MinOtf',1.2e-3,'lambdaEm',PupilInterpolators.lambda);
        psf{:}=VolumePSFSim(ImageParam,PSFParam);
    else
        error('Error generating PSF. Please supply the PSF parameters in the regularisation list for the PSF (3rd cell)');
    end
end
%% Normalize PSF
mysum=0;
if isreal(psf{1})
    ComplexPSF=0;
else
    ComplexPSF=1;
end
for v=1:numel(psf)
    if ~ComplexPSF  % complex valued PSFs are often 0 on average. It is thus normalized such that it's abssqr has a sum of one
        mysum=mysum+sum(psf{v});
    else
        mysum=mysum+sum(psf{v}*conj(psf{v}));
    end
end
for v=1:length(psf)
    if ~ComplexPSF  % complex valued PSFs are often 0 on average. It is thus normalized such that it's abssqr has a sum of one
        psf{v}=psf{v}/mysum;  % Forces the PSF to be normalized (over all sup-PSFs)
    else
        psf{v}=psf{v}/sqrt(mysum);  % Forces the PSF to be normalized 
    end
end
clear mysum;

otfrep=cell(numel(psf),1);


%%
for v=1:numel(psf)
    if ~isempty(psf{v})
    if ComplexPSF
        otfrep{v}=ft(psf{v});
    else
        % otfrep{v}=rft(dip_image(fftshift(double(psf{v})))); % The fft shift is necessary, since the rft assumes a different coordinate zero position
        otfrep{v}=rft(ifftshift(psf{v})); % The fft shift is necessary, since the rft assumes a different coordinate zero position
    end
    clear psf{v};  % These are not needed any longer. Free them, especially if converted to cuda
    if length(NewDataSize) > 2 && OrigSize(3) == 1 && (NewDataSize(3) == OrigSize(3))
       fprintf('Thick slice speedup: Adjusting OTF #%d for slicing operation.\n',v);
       if ~ ComplexPSF        % fft-based calculations are easier than rft calculations
           % the line below accounts for the fact for rft convolutions OTFs are centered at the border, not in the center, but a simle sum over Z is used
            otfrep{v}= otfrep{v}.* ifftshift(exp(i*(floor(size(otfrep{v},3)/2)/size(otfrep{v},3))*2*pi*zz(size(otfrep{v}))));
       end
       % otfrep{v}= otfrep{v} / sqrt(size(otfrep{v},3));    % This line turns out to be wrong. Nevertheless results can look nicer ??!
    end
    end
end

% if ~isempty(myillu_mask) %If there is a zero in the size of an array, cuda_cuda will give an error. Aurelie 03.03.2014
    myillu_mask=ConditionalCudaConvert(myillu_mask,useCuda,1);
% end
% if ~isempty(myillu) %If there is a zero in the size of an array, cuda_cuda will give an error. Aurelie 03.03.2014
    myillu=ConditionalCudaConvert(myillu,useCuda);
% end

mysum=0;
if ~isempty(DeconvMask)
    myDivisor=sum(DeconvMask,DeconvMask);
else
    myDivisor=prod(OrigSize);
end
measSums=cell(length(myim),1);
measSumsSqr=cell(length(myim),1);



for v=1:length(myim)
    measSums{v}=sum(myim{v},DeconvMask);
    if RegObj(17,1) || RegIllu(17,1) || RegOTF(17,1)
        measSumsSqr{v}=sqrt(sum(real(myim{v} .* conj(myim{v})),DeconvMask));
    end
    if ~isempty(myillu)
        mysum=mysum+measSums{v}/myDivisor/(mean(myillu{v}));  % use the old size because of the effect of zero padding and the normalisation of the psf
    else
        mysum=mysum+measSums{v}/myDivisor;  % use the old size because of the effect of zero padding and the normalisation of the psf
    end
end
mymean=mysum/length(myim);  % Force it to be real

    if (1)
        svecsize=NewIlluSize; % as defined above.  Not size(myim{1}) is the object lives in a different space
        if RegObj(7,1)  % this means a complex object shall be reconstructed. Thus the start vector also needs to be complex
            startVec=newim(svecsize,'scomplex');
        else
            startVec=newim(svecsize);
        end
        if RegObj(9,1) && (mymean < 1e-2)   % Force Positive. In this case one needs an offset, even though the model (e.g. in the amplite world) does not support it.
            startVec=startVec+1e-3;
            % startVec=NumViews*(double(repmat(1e-3,[1 prod(svecsize)])))';  % For now, just start with a guess of one, so the amplitude algorithm does not screw up
        else
            startVec=startVec+mymean-RegObj(12,1);  % account for the background value in the object estimate
            % startVec=NumViews*(double(repmat(mymean,[1 prod(svecsize)])))';
        end
    else
        global para;
        startVec=NumViews*para.res_object_resampled * mymean / mean(para.res_object_resampled)';   % To test it on SIM simulations
    end    
    
    if (~RegObj(7,1))  % Even though the data may be complex, the model may be forced to be real.
        startVec=abs(startVec); % sqrt will be applied later if needed. This is only to force it to be a real valued quantity even if the data is complex
    end
    if (RegObj(8,1)) % IntensityData
        for v=1:length(myim)
            if ~isreal(myim{1})
                error('For intensity data as specified in the flag, all images need to be real valued!');
            end
        end
        % startVec=startVec*exp(i*pi/4);  % Maybe random phases should be supplied as a start?
        if RegObj(7,1)  % 'complex'
            startVec=sqrt(startVec)+1e-8*i;  % Maybe random phases should be supplied as a start?
        else
            startVec=sqrt(startVec);  % Maybe random phases should be supplied as a start?
        end
        if ~RegObj(21,1)  % not FTData
            startVec=startVec/sqrt(sum(otfrep{1}));
        end
    end

AssignFunctions(RegularisationParameters,0); % Object estimate for the startVec estimation below.

if isempty(aRecon) || ~equalsizes(size(aRecon),size(startVec))
    aRecon=startVec;  % Just to have the size information inside
end

startVec=ConvertModelToVec(startVec);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
% startVec=convertGradToVec(startVec);  % Changes the format to Matlab and packs complex numbers appropriately

% Randomized starting vector helps sometime for Poisson, but seems bad as HF noise can end in Checkerboard patterns
% startVec=startVec + mymean*rand(size(startVec))/10;   % Introduce some extra noise. This sometimes helps the initial gradient

if ~isempty(DeconvMask)
        DeconvMask=ConditionalCudaConvert(DeconvMask,useCuda);
end
startVec=ConditionalCudaConvert(startVec,useCuda);

tic

% Advances speed up version:

% if lambda < 0  % choose lambda automatically
%     lambdaPenalty=0; % remove the penalty term for the line below
%     [err,grad]=GenericErrorAndDeriv(startVec);  % just to get an idea about the size of the gradient to expect
%     lambdaPenalty=7e4/(mean(mymean)/mean(abs(grad))); % 1e6;
%     fprintf('Lambda was estimated to %g\n',lambdaPenalty);
% end
clear mymean;

if ~isempty(RefObject)
        myScale = mean(RefObject) / mean(myim{1});
        if (abs(myScale -1.0) > 0.01)
            fprintf('WARNING! Reference Object and measured data have a significantly different mean value! Adapting the reference by a scaling factor %g\n',myScale);
            RefObject = RefObject / myScale;
        end
        if (~isvector(RefObject))
            RefObject=ConvertModelToVec(RefObject);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
        end
        % if ~RegObj(6,1)  % means reuse previous result
            RefObject_SSQ=[];  % For checking the progress to the ground truth during iterations.  
            RefObject_SAbs=[]; % For checking the progress to the ground truth during iterations.  
            RefObject_SSQ2=[];  % For checking the progress to the ground truth during iterations. 
            RefObject_SAbs2=[]; % For checking the progress to the ground truth during iterations. 
        % end
end

%lambdaPenalty=0; % 1e6;
if Update(1)~='B'   % not blind
    %NormFac=0.06;
    if (RegularisationParameters(18,1) <= 0)
            NormFac=1;            
            RegularisationParameters=RegObj;
            if isempty(myillu) || isempty(myillu{1})
                AssignFunctions(RegularisationParameters,0); % To estimate is Blind Object Step
            else
                AssignFunctions(RegularisationParameters,1); % To estimate is Blind Object Step
            end
            savedRecon=aRecon;
            [val,agrad]=GenericErrorAndDeriv(startVec);  % is used to determine a useful value of the normalisation
            aRecon=savedRecon; clear savedRecon;
            NormFac=1/(norm(agrad)/numel(agrad))/1e6; % /1e6
        %[val,agrad]=GenericErrorAndDeriv(startVec);  % is used to determine a useful value of the normalisation
        %aNorm=1/(norm(agrad)/numel(agrad));
        %NormFac=aNorm;
    else
        NormFac=RegularisationParameters(18,1); %0.06;  % 1e-6
    end
    fprintf('\nObject NormFac is %g\n',NormFac);
    if RegObj(6,1)  % means reuse previous result
        if ~equalsizes(size(aRecon),NewSize)
            fprintf('Warning! Start Data has wrong size. Adapting size.')
            aRecon=extract(aRecon,NewSize);
        end
        startVec=aRecon;
        startVec=ConvertModelToVec(startVec);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
    end
    RegularisationParameters=RegObj;
    if isempty(myillu)
        AssignFunctions(RegularisationParameters,0); % To estimate is Object (non-blind)
    else
        AssignFunctions(RegularisationParameters,1); % To estimate is Object but fixed illumination has to be considered
    end

    [myRes,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,startVec,NumIter(1));

    AssignToGlobal(ConvertInputToModel(myRes)); % Will save myRes to aRecon and correct for possible squaring ...
    if nargout > 3
          evolObj=myoutput.trace.fval / NormFac;
    end
else % 'B' Do a blind estimation of the unknown intensity or OTF distribution
%% Do some sanity checks here
if isempty(myillu_sumcond)
    myillu_sumcond={length(myim)};
end
if myillu_sumcond{1} < 0  % indicates to use no sum condition
    myillu_sumcond=[];
end
        
if ~isempty(myillu_sumcond)
    for si= 1:size(myillu_sumcond,1)  % sub illuminations. Aurelie 23.06.2014
    for n=1:size(myillu_sumcond,2) % Aurelie 23.06.2014. Replace myillu_sumcond{n} by myillu_sumcond{si,n} in this loop
        if myillu_sumcond{si,n}<2 && (length(myim) > 1)
            error('First sum condition has to be at least two.');
        end
%         if myillu_sumcond{n}>numel(myim)
        if myillu_sumcond{si,n}>numel(myim) && (numel(otfrep)==1) %Aurelie 27.05.2014: in the case of multi-view deconvolution, the sum condition may be larger
            error('Last sum condition cannot be larger than number of images.');
        end
        if n>1 && (myillu_sumcond{si,n} - myillu_sumcond{si,n-1}) < 2
            error('Distance between sum conditions has to be at least two.');
        end
    end
    end
end
        Update=Update(2:end);  % remove the first letter (standing for "Blind")
        % Update='lbfgs'  % 'cg'
        DataSize=prod(size(myim{1}));

        if RegIllu(6,1)
            VecIllu=myillu;
        elseif NumIter(4) ~= 0 % ~isempty(Regularisation{2})  % This means that a spatially variing illumination is part of the model. Thus a start vector is needed
            if useCuda
%                 VecIllu = repmat(cuda(dip_image(1)),[NewIlluSize numel(myim)-length(myillu_sumcond)]);
%                 VecIllu = repmat(cuda(dip_image(1)),[NewIlluSize (numel(myim)-length(myillu_sumcond))*numel(otfrep)]); %Aurelie 27.05.2014
                VecIllu = repmat(cuda(dip_image(1)),[NewIlluSize (numel(myim)-length(myillu_sumcond))]); %Aurelie 23.06.2014
            else
                VecIllu = repmat(dip_image(1),[NewIlluSize (numel(myim)-length(myillu_sumcond))*numel(otfrep)]); %Aurelie 27.05.2014
            end
            %         if isa(myim{1},'cuda')
            %             VecIllu=ones_cuda(DataSize*(numel(myim)-length(myillu_sumcond)),1);  % one less than all illumination patterns, as the last one is defined by the sum condition
            %         else
            %             VecIllu=ones(DataSize*(numel(myim)-length(myillu_sumcond)),1);  % one less than all illumination patterns, as the last one is defined by the sum condition
            %         end
            %
            AssignFunctions(RegularisationParameters,2); % Illumination estimate for the startVec estimation below.
            VecIllu=ConvertModelToVec(VecIllu);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
            %
            sumSqr=0;
            for v=1:length(myim)
                myillu{v}=newim(NewIlluSize)+1;  % start with uniform illumination
                %myillu{v}=myim{v}*0+2*rand(size(myim{1},1),size(myim{1},1));  % start with uniform illumination
                sumSqr=sumSqr+sum(myim{v}.*myim{v});   % for calculating the normalisation
                % VecIllu(1+DataSize*(v-1):DataSize*v)=(double(reshape(myillu{v},[DataSize 1])))';
            end
        end
        %% estimate a good NormFac to use
        if (RegularisationParameters(18,1) <= 0)
            NormFac=1;
            RegularisationParameters=RegObj;AssignFunctions(RegularisationParameters,1); % To estimate is Blind Object Step
            savedRecon=aRecon;
            [msevalue,agrad]=GenericErrorAndDeriv(startVec);  % is used to determine a useful value of the normalisation
            aRecon=savedRecon; clear savedRecon;
            aNorm=1/(norm(agrad)/numel(agrad))  / 100;
        else
            aNorm=RegularisationParameters(18,1); % 0.06;
            % aNorm=0.0001;
        end
        NormFac=aNorm;
        fprintf('Object NormFac is %g\n',NormFac);
        if RegObj(6,1)  % means reuse previous result
            startVec=aRecon;
            startVec=ConvertModelToVec(startVec);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
        end
        NormFacOld=aNorm;
        
        myRes=startVec;
        AssignFunctions(RegularisationParameters,0); % Object estimate for the startVec estimation below.
        AssignToGlobal(ConvertInputToModel(myRes)); % ,size(myim{1})  % will square the startvec, (if needed) to produce aRecon

        evolObj=[];
        evolIllu=[];
        NumObjIter=0; 
        NumIlluIter=NumIter(4); 

        for n=1:NumIter(1)
            ToEstimate=0;  % object estimation step
            %NormFac=0.06;
            NormFac=aNorm;
            fprintf('Object NormFac is %g, Prev. Error was %g\n',NormFac,msevalue*NormFac/NormFacOld);
            NormFacOld=NormFac;
            %NormFac=1/sumSqr;
            InstantUpdate=1;  % defines whether the object update is done directly after the object iteration or only after illumination has been updated as well.
            if n==1
                NumObjIter=NumIter(2); % SIM 5;  Speckles: 0
            else
                NumObjIter=NumIter(3); %SIM 25;  Speckles: 5
            end

            savedRecon=aRecon; % make sure it does not get overwritten
            fprintf('\n\nIterating Object, cycle %d\n',n);
            if NumObjIter>0
                RegularisationParameters=RegObj;
                AssignFunctions(RegularisationParameters,1); % To estimate is Blind Object Step
                [myRes,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,myRes,NumObjIter);  % Will also alter aRecon
                if nargout > 3
                    evolObj =[evolObj ; myoutput.trace.fval/ NormFac];
                end
                if InstantUpdate
                    AssignToGlobal(ConvertInputToModel(myRes)); % Will save myRes to aRecon
                    dipshow(3,aRecon);drawnow();
                else
                    aRecon=savedRecon; % make sure it does not get overwritten
                end
            end
            %[myRes,msevalue,moreinfo]=minFunc(@GenericErrorAndDeriv,myRes,optionsObj); % @ means: 'Function handle creation' 
            % aRecon should still be globally set

            %saveRecon=aRecon;
            if NumIlluIter ~= 0
                ToEstimate=1;  % illumination estimation step
                % aRecon=TheObject;
                if n==1
                    RegularisationParameters=RegIllu;AssignFunctions(RegularisationParameters,2); % To estimate is Blind Illumination Step
                    if (RegularisationParameters(18,1) <= 0) % can also be if (1)
                        NormFac=1;
                        [val,agrad]=GenericErrorAndDeriv(VecIllu);  % is used to determine a useful value of the normalisation
                        IlluNorm=1/(norm(abs(agrad))/numel(agrad))/1000;
                    else
                        IlluNorm=RegularisationParameters(18,1);
                    end
                    %[val,agrad]=GenericErrorAndDeriv(VecIllu);  % is used to determine a useful value of the normalisation
                    %(norm(agrad)/numel(agrad))
                end
                fprintf('Illu NormFac is %g, Prev. Error was %g\n',IlluNorm,msevalue*IlluNorm/NormFacOld);
                NormFacOld=IlluNorm;
                %NormFac=IlluNorm;
                %% Illumination estimation step
                if size(NumIter,2) > 3
                    % NormFac=1/sumSqr;
                    % NormFac=0.1; % /sumSqr;
                    fprintf('\n\nIterating Illumination, cycle %d\n',n);
                    RegularisationParameters=RegIllu;AssignFunctions(RegularisationParameters,2); % To estimate is Blind Illumination Step
%                     NormFac=IlluNorm*6e-5;
                    NormFac=IlluNorm; 
                    [VecIllu,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,VecIllu,NumIlluIter);
                    if nargout > 4 && ~isempty(myoutput)
                        evolIllu=[evolIllu ; myoutput.trace.fval / NormFac];
                    end
                    %[VecIllu,msevalue,moreinfo]=minFunc(@GenericErrorAndDeriv,VecIllu,optionsIllu); % @ means: 'Function handle creation'
                    AssignToGlobal(ConvertInputToModel(VecIllu));  % writes result into the myillu images
                    illu=cat(4,myillu{:});
                    dipshow(4,illu);drawnow();
                    if nargout > 7
                        allIllu{n}=illu;
                    end
                    clear illu;
                    if ~InstantUpdate
                        AssignFunctions(RegularisationParameters,1); % To estimate is Blind Illumination Step
                        AssignToGlobal(ConvertInputToModel(myRes));
                        dipshow(3,aRecon);drawnow();
                    end
                end
            end
            %% Estimating the OTF
            ToEstimate=2;  % OTF estimation step
            if size(NumIter,2) > 4  && NumIter(5) ~= 0
                % NormFac=1e-10/sum(aRecon);
                NumOTFIter=NumIter(5); % SIM 5;  Speckles: 15
                fprintf('\n\nIterating PSF, cycle %d\n',n);
                RegularisationParameters=RegOTF;
                AssignFunctions(RegularisationParameters,3); % To estimate is Blind Illumination Step (OTF estimation)
                if n==1
                    if (1) % isempty(OTFmask) || isempty(OTFmask{v}) % creat OTF masks from ideal (unaberrated) OTF if necessary
                        OTFmask=cell(numel(otfrep),1);
                        for no=1:numel(otfrep)
                            absotf=abs(otfrep{no});
                        	OTFmask{no}=absotf>(max(absotf)/10000);
                        end
                    end
                    allotf=cat(4,otfrep{:});
                    % flag below is noFFT and allows an OTF instead of PSF to be used here.
                    global noFFT;  % This is a hack to prevent the fft as this is already an OTF
                    noFFT=1;
                    if RegularisationParameters(14,1)  % This means the 2D pupil ist really what is estimated
                        mysize=size(PupilInterpolators.indexList2D,2);
                        VecOTF=ones(2*mysize,1)/sqrt(2);  %/mysize
                        ResOTF=ConvertInputToModel(VecOTF);
                        VecOTF=VecOTF/sqrt(sum(rift(ResOTF{1})));  % To force the PSF to have an integral of one
%                        VecOTF=double(VecOTF(:));
                    else  % estimate the 3D OTF voxels
                        VecOTF=ConvertModelToVec(allotf);    % converts the OTF (not PSF!) dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
                        if iscell(VecOTF)
                            VecOTF=VecOTF{1};
                        end
                    end
                    startVec=VecOTF;
                    noFFT=0;
                    if (RegularisationParameters(18,1) <= 0) % can also be if (1)
                        NormFac=1;
                        [val,agrad]=GenericErrorAndDeriv(startVec);  % is used to determine a useful value of the normalisation
                        NormFacOTF=1/(norm(agrad)/numel(agrad))*1e-8;
                    else
                        NormFacOTF=RegularisationParameters(18,1);
                    end
                end
                NormFac=NormFacOTF;
                fprintf('OTF NormFac is %g, Prev. Error was %g\n',NormFacOTF,msevalue*NormFac/NormFacOld);
                NormFacOld=NormFac;
                [VecOTF,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,VecOTF,NumOTFIter);
                if nargout > 3
                    evolOTF=[evolOTF ; myoutput.trace.fval];
                end
                %[VecIllu,msevalue,moreinfo]=minFunc(@GenericErrorAndDeriv,VecIllu,optionsIllu); % @ means: 'Function handle creation'
                AssignToGlobal(ConvertInputToModel(VecOTF));  % writes result into the otfrep images
                if nargout > 6
                    allOTF{n}=otfrep;
                end
                allotf=cat(4,otfrep{:});
                psf=fftshift(rift(otfrep{1}));
                dipshow(6,SubSlice(psf,2)); % (:,floor(size(psf,2)/2),:)
                global savedATF;
                if ~isempty(savedATF)
                    dipshow(7,savedATF);
                    dipshow(8,phase(savedATF));
                end
                dipshow(5,allotf);drawnow();
            end

        end
    if (NumObjIter > 0)    
        AssignFunctions(RegObj,0); % To estimate is Blind Illumination Step
        AssignToGlobal(ConvertInputToModel(myRes));
    end
    if (NumIlluIter > 0)    
        AssignFunctions(RegIllu,2); % To estimate is Blind Illumination Step
        AssignToGlobal(ConvertInputToModel(VecIllu));
    end
    if size(NumIter,2) > 4  && NumIter(5) ~= 0 %Aurelie 25.03.2014: NumOTFIter is not assigned   
        AssignFunctions(RegOTF,3); % To estimate is Blind Illumination Step
        AssignToGlobal(ConvertInputToModel(VecOTF));
    end
end
mytime=toc;
%fprintf('Time was %g, lambda: %g\n', mytime, lambdaPenalty)
fprintf('Time was %g\n', mytime)

% res=reshape(dip_image(myRes','single'),size(myim{1}));
%if isa(myRes,'cuda')
%    convertVecToObj(double_force(myRes));
%els
%end

clear myim;
clear DeconvMask;
if useCuda
    cuda_clearheap();
    toc
    aRecon=dip_image_force(aRecon);
    toc
    set_ones_cuda(0); set_zeros_cuda(0);   % make sure that the zeros and ones function from matlab is used.
end
res=aRecon;
clear aRecon;


if (nargout > 1)
    resIllu=myillu;
end

if (nargout > 2)
    resPSF=cell(1,numel(otfrep));
    for c=1:numel(otfrep)
        resPSF{c}=ifftshift(rift(otfrep{c}));
    end
end
clear otfrep;

if norm(borderSizes) > 0 && keepExtendedStack == 0
   res=extract(res,subSampling .* floor(aResampling .* OrigSize(1:length(size(res)))));
   fprintf('Reduced size to original\n');
   if  exist('resPSF')
       for d=1:length(resPSF)
           resPSF{d}=extract(resPSF{d},floor(aResampling .* OrigSize(1:length(size(res)))));
       end
   end
end

if nargout > 6
    allRObj = allObj;
end

global FinalErrNorm
FinalErrNorm=msevalue/NormFac;
fprintf('Normalized Error norm is %g\n',FinalErrNorm)

if ~isempty(RefObject)
    figure
    plot(RefObject_SAbs,'b');hold on
    plot(RefObject_SSQ,'g');
    title('Estimation vector')
    legend({'Mean Abs Error','Sqrt(Mean Squared Error)'})
    ylabel('Error to ground truth');
    xlabel('Iteration number');    
    fprintf('To access the error values type\nglobal RefObject_SAbs;global RefObject_SSQ;\n')
    figure
    plot(RefObject_SAbs2,'b');hold on
    plot(RefObject_SSQ2,'g');
    title('Squared estimation vector')
    legend({'Mean Abs Error','Sqrt(Mean Squared Error)'})
    ylabel('Error to ground truth');
    xlabel('Iteration number');    
    fprintf('To access the error values type\nglobal RefObject_SAbs2;global RefObject_SSQ2;\n')

    figure
    plot(RefObject_SAbs/RefObject_SAbs(1),'b');hold on
    plot(RefObject_SSQ/RefObject_SSQ(1),'g');
    title('Estimation vector')
    legend({'Mean Abs Error','Sqrt(Mean Squared Error)'})
    ylabel('Normalized Error to ground truth');
    xlabel('Iteration number');    
    figure
    plot(RefObject_SAbs2/RefObject_SAbs2(1),'b');hold on
    plot(RefObject_SSQ2/RefObject_SSQ2(1),'g');
    title('Squared estimation vector')
    legend({'Mean Abs Error','Sqrt(Mean Squared Error)'})
    ylabel('Normalized Error to ground truth');
    xlabel('Iteration number');    

    RefObject=[];  % To not always track stuff
end
