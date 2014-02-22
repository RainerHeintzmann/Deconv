% [res,resIllu,resPSF,evolObj,evolIllu,evolOTF]=GenericDeconvolution(image,psf,NumIter,Method,Update,Regularisation,betas,borderSizes,Variances,useCuda,keepExtendedStack) : Deconvolution routine with various functionalities.
% image: Image to Deconvolve
% psf: Point spread function to deconvolve with
% NumIter: Iterations to perform. Optionally the subiterations for initial Object iterations (2) and successive object (3), illumination (4) and psf (5) iterations  can be given.
% Method: One of 'LeastSqr', 'WeightedLeastSqr' and 'Poisson', defining the norm of the data-simulation agreement to optimize
% Update: Default: [] uses 'lbfgs', This describes the minimization algorithm (update scheme) to use. You can use all the ones that minFunc allows, but in addition also 'RL' for
% RichardsonLucy multiplicative (EM) Algorithm, or 'RLL' for a Wolfe line search along the Richardson Lucy update direction
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
%       'AR',[lambda,gamma]: EM-regularization according to Arigovindan et al., PNAS, 2013. penalty=ln(epsR+ (|f|^2+ SumHessianSqr)),  SumHessianSqr= (d^2/dxdx f)^2 + (d^2/dydy f)^2+ (2d^2/dxdy f)^2
%       'Bg', value: Fixed background value (to get the Poisson Statistics right)
%       'ForcePos',[] : Will make the object the square of an auxiliary function, which is then estimated
%       'Complex',[] : Will allow the object to be an aplitude object. Note that even for complex valued data the object can remain to be real.
%       'Resample',[sx sy sz] : Specifies the resampling to be used for reconstructed object (only for object at the moment)
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
% Example calls:
% For Poisson noise, Goods roughness and positivity constraint and twice resampling of the output:
%     myDeconvGP=GenericDeconvolution(img,h,85,'Poisson',[],{'GR',0.01;'ForcePos',[];'Resample',2},[1,1,1],[0 0 0],[],useCuda); gtp=cat(1,img{1},myDeconvGP)
% For complex data and complex PSF
%     myDeconv=GenericDeconvolution(img,h,15,'LeastSqr',[],{'Complex',[]},[1 1 1],[0 0 0],[],useCuda); st=cat(1,img{1},myDeconv)


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

function [res,resIllu,resPSF,evolObj,evolIllu,evolOTF]=GenericDeconvolution(image,psf,NumIter,Method,Update,Regularisation,betas,borderSizes,Variances,useCuda,keepExtendedStack)
global myim; % Cell array of the images
global otfrep; % Cell array of the images
global DeconvMethod;
global aResampling; 
global RegularisationParameters;
global DeconvMask;
global DeconvVariances;
global BetaVals;
global myillu; 
global myillu_mask; 
global OTFmask;
global ToEstimate;   % This flag controls what to estimate in this iteration step. 0: sample density, 1: illumination intensity, 2: both, 3: psf
global aRecon;   % This flag controls what to estimate in this iteration step. 0: sample density, 1: illumination intensity, 2: both, 3: psf
global NormFac;  % normalisation factor
global myillu_sumcond;
% global TheObject;
global cuda_enabled;  % also to see if cuda is installed
global useCudaGlobal;useCudaGlobal=useCuda;
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF
global ConvertModelToVec; % Contains the routine to convert the model (e.g. the gradient of the object) into the vector to iterate

%global ComplexObj;    % Can the object to estimate be complex valued?
%global IntensityData;  % Is the supplied data the abs square of the object to estimate?
%global ForcePos;
%ForcePos=0;         % Will only be set to one in the ParseRegularisation routine, if chosen.
global ComplexPSF; % will be set, if PSF is complex valued
DeconvMethod=Method;
aResampling=1; 

if nargin < 11
    keepExtendedStack=0;
end
if nargin < 10
    useCuda=0;
end
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
[RegObj,RegIllu,RegOTF]=ParseRegularisation(Regularisation);  % Also sets global variable such as ForcePos and ComplexObj
RegularisationParameters=RegObj;
% if isempty(myillu)
%     AssignFunctions(RegularisationParameters,0); % To estimate is Object (non-blind)
% else
%     AssignFunctions(RegularisationParameters,1); % To estimate is Object but fixed illumination has to be considered
% end

%%
% somehow the line below gets overwritten, so we have to set the boudary
% in the inner loop "GenericErrorAndDeriv"
dip_setboundary('periodic')  % Establishes Periodic Boudary conditions.

borderSizes=borderSizes*2;

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
        
    extraBorder(2) = mod(s(2)+borderSizes(2),2);  % Only expand for uneven sizes along Y
    if numel(borderSizes)> 2 && borderSizes(3)==0
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
    if useCuda
        for v=1:length(myim)
            myim{v}=cuda(myim{v});
        end
    end
else
    for v=1:size(image,4)
        if length(size(image))>3
            myim{v}=image(:,:,:,v-1);
        else
            myim{v}=image;
        end
        if useCuda
            myim{v}=cuda(myim{v});
        end
    end
end
clear image;

if ~(iscell(psf) || isa(psf,'dip_image_array'))
    for v=1:size(psf,4)
        if length(size(psf))>3
            n{v}=psf(:,:,:,v-1);
        else
            n{v}=psf;
        end
        if useCuda
            n{v}=cuda(n{v});
        end
    end
    psf=n;
else
   if useCuda
        for v=1:length(psf)
            psf{v}=cuda(psf{v});
        end
    end

end

OrigSize=size(myim{1});
if length(borderSizes) > length(OrigSize)
   OrigSize(end+1:length(borderSizes))=1;
end
%% Extend the boarder
if isa(borderSizes,'dip_image')
    fprintf('Mask provided as bordersize. Using mask instead of generating borders');
    DeconvMask=borderSizes;    
else if norm(borderSizes) > 0
        NewSize=floor(OrigSize+borderSizes);
        fprintf('Border region requested, expanding size from [');
        fprintf('%g,',OrigSize);
        fprintf('] to [');
        fprintf('%g,',NewSize);
        fprintf('] \n');
        DeconvMask=extract(myim{1}*0>-1,NewSize);  % Place zero outside and one at the old data region
        if ~isempty(Variances)
            DeconvVariances=extract(Variances,NewSize,floor(NewSize/2),1);  % Pad with Variance of one (gets ignored)
        end
        for v=1:length(myim)
            myim{v}=extract(myim{v},NewSize);  % Overwrite the old data            
            if (v<=numel(psf)) && norm(size(psf{v})-NewSize)~=0
               psf{v}=squeeze(extract(psf{v},NewSize,[],'cyclic'));
               fprintf('Warning: View %d, PSF has incorrect size. Adapting the size of the PSF\n',v);
            end
            if exist('myillu') && ~isempty(myillu)
                sillu=size(myillu{1+mod(v-1,numel(myillu))});
                sillu(end+1:length(NewSize))=1;
                if (v <= numel(myillu)) && norm(sillu-NewSize) ~= 0
                    fprintf('Warning: View %d, The given illumination distribution does not correspond in size to the extended border region. Cyclically expanding to fix this.\n',v)
                    myillu{v}=extract(myillu{v},NewSize,[],'cyclic');
                end
            end
        end
    else
        for v=1:numel(psf)
            spsf=size(psf{1+mod(v-1,numel(psf))});
            OrigSize(end+1:length(spsf))=1;
            if norm(spsf-OrigSize)~=0
                psf{1+mod(v-1,numel(psf))}=squeeze(extract(psf{v},OrigSize));   % In case the PSF has a different size
               fprintf('Warning: View %d, PSF has incorrect size. Adapting the size of the PSF\n',v);
            end
        end
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
    if useCuda
       if ComplexPSF
           otfrep{v}=ft(cuda(psf{v}));  
       else
           otfrep{v}=rft(fftshift(cuda(psf{v})));  % The fft shift is necessary, since the rft assumes a different coordinate zero position
       end
    else
       if ComplexPSF
           otfrep{v}=ft(psf{v});
       else
           % otfrep{v}=rft(dip_image(fftshift(double(psf{v})))); % The fft shift is necessary, since the rft assumes a different coordinate zero position
           otfrep{v}=rft(fftshift(psf{v})); % The fft shift is necessary, since the rft assumes a different coordinate zero position
       end
    end
    clear psf{v};  % These are not needed any longer. Free them, especially if converted to cuda
end

if useCuda
    for v=1:numel(myillu_mask)
        if ~isa(myillu_mask{v},'cuda')
            myillu_mask{v}=cuda(myillu_mask{v});
        end
    end
    for v=1:numel(myillu)
        if ~isa(myillu{v},'cuda')
            myillu{v}=cuda(myillu{v});
        end
    end
else
    for v=1:numel(myillu_mask)
        if isa(myillu_mask{v},'cuda')
            myillu_mask{v}=dip_image_force(myillu_mask{v}) ~= 0;
        end
    end
    for v=1:numel(myillu)
        if isa(myillu{v},'cuda')
            myillu{v}=dip_image_force(myillu{v});
        end
    end
end

mysum=0;
if ~isempty(DeconvMask)
    myDivisor=sum(DeconvMask,DeconvMask);
else
    myDivisor=prod(OrigSize);
end
for v=1:length(myim)
    if ~isempty(myillu)
        mysum=mysum+sum(myim{v},DeconvMask)/myDivisor/(sum(myillu{v})/prod(OrigSize));  % use the old size because of the effect of zero padding and the normalisation of the psf
    else
        mysum=mysum+sum(myim{v},DeconvMask)/myDivisor;  % use the old size because of the effect of zero padding and the normalisation of the psf
    end
end
mymean=mysum/length(myim);  % Force it to be real

if RegObj(6,1)  % means reuse previous result
    startVec=aRecon;
else
    if (1)
        svecsize=size(myim{1});
        if any(aResampling~=1)
            svecsize = floor(svecsize .* aResampling);
        end
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
        aRecon=startVec;
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
        if RegObj(7,1)
            startVec=sqrt(startVec)+1e-8*i;  % Maybe random phases should be supplied as a start?
        else
            startVec=sqrt(startVec);  % Maybe random phases should be supplied as a start?
        end
        startVec=startVec/sqrt(sum(otfrep{1}));
    end
end

AssignFunctions(RegularisationParameters,0); % Object estimate for the startVec estimation below.
startVec=ConvertModelToVec(startVec);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
% startVec=convertGradToVec(startVec);  % Changes the format to Matlab and packs complex numbers appropriately

% Randomized starting vector helps sometime for Poisson, but seems bad as HF noise can end in Checkerboard patterns
% startVec=startVec + mymean*rand(size(startVec))/10;   % Introduce some extra noise. This sometimes helps the initial gradient

if useCuda
    if ~isempty(DeconvMask)
        DeconvMask=cuda(DeconvMask);
    end
    startVec=cuda(startVec);
    set_ones_cuda(useCuda); set_zeros_cuda(useCuda);   % make sure that the zeros and ones function from cuda is used.
end

tic

% Advances speed up version:

% if lambda < 0  % choose lambda automatically
%     lambdaPenalty=0; % remove the penalty term for the line below
%     [err,grad]=GenericErrorAndDeriv(startVec);  % just to get an idea about the size of the gradient to expect
%     lambdaPenalty=7e4/(mean(mymean)/mean(abs(grad))); % 1e6;
%     fprintf('Lambda was estimated to %g\n',lambdaPenalty);
% end
clear mymean;

%lambdaPenalty=0; % 1e6;
if Update(1)~='B'
     NormFac=0.06;
%     [val,agrad]=GenericErrorAndDeriv(startVec);  % is used to determine a useful value of the normalisation
%     aNorm=1/(norm(agrad)/numel(agrad));
%     NormFac=aNorm;
    % NormFac=0.06;  % 1e-6

    RegularisationParameters=RegObj;
    if isempty(myillu)
        AssignFunctions(RegularisationParameters,0); % To estimate is Object (non-blind)
    else
        AssignFunctions(RegularisationParameters,1); % To estimate is Object but fixed illumination has to be considered
    end

    [myRes,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,startVec,NumIter(1));

    ConvertInputToModel(myRes); % Will save myRes to aRecon and correct for possible squaring ...
    if nargout > 2
          evolObj=myoutput.trace.fval;
    end
else % Do a blind estimation of the unknown intensity or OTF distribution
%% Do some sanity checks here
if ~isempty(myillu_sumcond)
    for n=1:numel(myillu_sumcond)
        if myillu_sumcond{n}<2 && (length(myim) > 1)
            error('First sum condition has to be at least two.');
        end
        if myillu_sumcond{n}>numel(myim)
            error('Last sum condition cannot be larger than number of images.');
        end
        if n>1 && (myillu_sumcond{n} - myillu_sumcond{n-1}) < 2
            error('Distance between sum conditions has to be at least two.');
        end
    end
end
        Update=Update(2:end);  % remove the first letter (standing for "Blind")
        % Update='lbfgs'  % 'cg'
        DataSize=prod(size(myim{1}));
        if isempty(myillu_sumcond)
            myillu_sumcond={length(myim)};
        end

        if RegIllu(6,1)
            VecIllu=myillu;
        elseif NumIter(4) > 0 % ~isempty(Regularisation{2})  % This means that a spatially variing illumination is part of the model. Thus a start vector is needed
            VecIllu = repmat(myim{1}*0+1,[1 1 1 numel(myim)-length(myillu_sumcond)]);
            %         if isa(myim{1},'cuda')
            %             VecIllu=ones_cuda(DataSize*(numel(myim)-length(myillu_sumcond)),1);  % one less than all illumination patterns, as the last one is defined by the sum condition
            %         else
            %             VecIllu=ones(DataSize*(numel(myim)-length(myillu_sumcond)),1);  % one less than all illumination patterns, as the last one is defined by the sum condition
            %         end
            %
            AssignFunctions(RegularisationParameters,2); % Illumination estimate for the startVec estimation below.
            VecIllu=ConvertModelToVec(VecIllu);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
            AssignFunctions(RegularisationParameters,0); % Illumination estimate for the startVec estimation below.
            %
            sumSqr=0;
            for v=1:length(myim)
                myillu{v}=myim{v}*0+1;  % start with uniform illumination
                %myillu{v}=myim{v}*0+2*rand(size(myim{1},1),size(myim{1},1));  % start with uniform illumination
                sumSqr=sumSqr+sum(myim{v}.*myim{v});   % for calculating the normalisation
                % VecIllu(1+DataSize*(v-1):DataSize*v)=(double(reshape(myillu{v},[DataSize 1])))';
            end
        end
        %% estimate a good NormFac to use
        if (1)
            NormFac=1;
            RegularisationParameters=RegObj;AssignFunctions(RegularisationParameters,1); % To estimate is Blind Object Step
            [val,agrad]=GenericErrorAndDeriv(startVec);  % is used to determine a useful value of the normalisation
            aNorm=1/(norm(agrad)/numel(agrad));
            fprintf('Object NormFac is %g\n',aNorm);
        else
            aNorm=0.06;
            % aNorm=0.0001;
        end
        NormFac=aNorm;
        
        myRes=startVec;
        ConvertInputToModel(myRes); % ,size(myim{1})  % will square the startvec, (if needed) to produce aRecon

        evolObj=[];
        evolIllu=[];
        for n=1:NumIter(1)
            ToEstimate=0;  % object estimation step
            %NormFac=0.06;
            NormFac=aNorm;
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
                if nargout > 2
                    evolObj =[evolObj ; myoutput.trace.fval];
                end
            end
            %[myRes,msevalue,moreinfo]=minFunc(@GenericErrorAndDeriv,myRes,optionsObj); % @ means: 'Function handle creation' 
            % aRecon should still be globally set
            if InstantUpdate
                ConvertInputToModel(myRes); % Will save myRes to aRecon
                dipshow(3,aRecon);drawnow();
            else
                aRecon=savedRecon; % make sure it does not get overwritten
            end
            %saveRecon=aRecon;
            if NumIter(4) > 0
                ToEstimate=1;  % illumination estimation step
                NumIlluIter=NumIter(4); % SIM 5;  Speckles: 15
                % aRecon=TheObject;
                if n==1
                    NormFac=1;
                    RegularisationParameters=RegIllu;AssignFunctions(RegularisationParameters,2); % To estimate is Blind Illumination Step
                    [val,agrad]=GenericErrorAndDeriv(VecIllu);  % is used to determine a useful value of the normalisation
                    IlluNorm=1/(norm(abs(agrad))/numel(agrad));
                    fprintf('Illu NormFac is %g\n',IlluNorm);
                    %[val,agrad]=GenericErrorAndDeriv(VecIllu);  % is used to determine a useful value of the normalisation
                    %(norm(agrad)/numel(agrad))
                end
                %NormFac=IlluNorm;
                %% Illumination estimation step
                if size(NumIter,2) > 3
                    % NormFac=1/sumSqr;
                    % NormFac=0.1; % /sumSqr;
                    fprintf('\n\nIterating Illumination, cycle %d\n',n);
                    RegularisationParameters=RegIllu;AssignFunctions(RegularisationParameters,2); % To estimate is Blind Illumination Step
                    NormFac=IlluNorm*6e-5;
                    [VecIllu,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,VecIllu,NumIlluIter);
                    if nargout > 2
                        evolIllu=[evolIllu ; myoutput.trace.fval];
                    end
                    %[VecIllu,msevalue,moreinfo]=minFunc(@GenericErrorAndDeriv,VecIllu,optionsIllu); % @ means: 'Function handle creation'
                    ConvertInputToModel(VecIllu);  % writes result into the myillu images
                    illu=cat(4,myillu{:});
                    dipshow(4,illu);drawnow();
                    clear illu;
                    if ~InstantUpdate
                        AssignFunctions(RegularisationParameters,1); % To estimate is Blind Illumination Step
                        ConvertInputToModel(myRes);
                        dipshow(3,aRecon);drawnow();
                    end
                end
            end
            %% Estimating the OTF
            ToEstimate=2;  % OTF estimation step
            if size(NumIter,2) > 4  && NumIter(5) > 0
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
                        global PupilInterpolators;
                        mysize=size(PupilInterpolators.indexList2D,2);
                        VecOTF=ones(2*mysize,1)/sqrt(2);  %/mysize
                        VecOTF=VecOTF/sqrt(sum(rift(ConvertInputToModel(VecOTF))));  % To force the PSF to have an in
                    else  % estimate the 3D OTF voxels
                        VecOTF=ConvertModelToVec(allotf);    % converts the OTF (not PSF!) dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
                    end
                    startVec=VecOTF;
                    noFFT=0;
                    NormFac=1;
                    [val,agrad]=GenericErrorAndDeriv(startVec);  % is used to determine a useful value of the normalisation
                    NormFacOTF=1/(norm(agrad)/numel(agrad));
                    fprintf('OTF NormFac is %g\n',NormFacOTF);
                end
                NormFac=NormFacOTF*1e-8;
                [VecOTF,msevalue,moreinfo,myoutput]=DoDeconvIterations(Update,VecOTF,NumOTFIter);
                if nargout > 3
                    evolOTF=[evolOTF ; myoutput.trace.fval];
                end
                %[VecIllu,msevalue,moreinfo]=minFunc(@GenericErrorAndDeriv,VecIllu,optionsIllu); % @ means: 'Function handle creation'
                ConvertInputToModel(VecOTF);  % writes result into the otfrep images
                allotf=cat(4,otfrep{:});
                dipshow(5,allotf);drawnow();
            end

        end            
end
%set_ones_cuda(0);
%set_zeros_cuda(0);
mytime=toc;
%fprintf('Time was %g, lambda: %g\n', mytime, lambdaPenalty)
fprintf('Time was %g\n', mytime)

% res=reshape(dip_image(myRes','single'),size(myim{1}));
%if isa(myRes,'cuda')
%    convertVecToObj(double_force(myRes));
%els
    AssignFunctions(RegularisationParameters,1); % To estimate is Blind Illumination Step
    ConvertInputToModel(myRes);
%end

clear myim;
clear DeconvMask;
if useCuda
    cuda_clearheap();
    toc
    aRecon=dip_image_force(aRecon);
    toc
    set_ones_cuda(0); set_zeros_cuda(0);   % make sure that the zeros and ones function from cuda is used.
end
res=aRecon;
clear aRecon;

if norm(borderSizes) > 0 && keepExtendedStack == 0
   res=extract(res,floor(aResampling .* OrigSize(1:length(size(res)))));
   fprintf('Reduced size to original\n');
end

if (nargout > 1)
    resIllu=myillu;
end

if (nargout > 1)
    resPSF=cell(1,numel(otfrep));
    for c=1:numel(otfrep)
        resPSF{c}=ifftshift(rift(otfrep{c}));
    end
end
clear otfrep;

