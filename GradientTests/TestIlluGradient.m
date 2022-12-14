%03.03.2014: Test new deconv toolbox OK. Aurelie

% This function simulates a small image and test numerically whether the gradiests as estimated by GenericErrorAndDrivative correspond to what is expected from the error value return of that function
%clear all
%disableCuda();
useCuda=1; disableCuda();
useSumCond=1;
sX=10; %10 and 50: checked
sY=10;
NumIm=2;
if (0)
    MyReg={};  % Why does this not work?
else
    MyReg={'ForcePos',[]};  % Why does this not work?
end

rng(1); % initialize the raqndom generator with the same seed always
obj=dip_image(rand(sX,sY));

%%estimate=mean(img)+img*0;  % Careful: This yields rubbish with estimating Regularisations
%%estimate=obj*1.1;
%%estimate=newim(9,9);estimate(4,4)=0.5;

%obj=newim(9,9);obj(4,4)=1;
% estimate=dip_image(rand(sX,sY)+i*rand(sX,sY));  % new OTF estimate (everywhere, not only in the mask)
rng(2); % initialize the raqndom generator with the same seed always
objestimate=dip_image(rand(sX,sY));  % current (old) reconstruction estimate
rng(3); % initialize the raqndom generator with the same seed always
illu=dip_image(rand(sX,sY,NumIm));   % an illumination distribution
rng(4); % initialize the raqndom generator with the same seed always
estimate=dip_image(rand(sX,sY));
if (0) %calculate with the mask
    global myillu_mask;
    myillu_mask=cell(1,NumIm);
    tmp=newim(estimate);
    mykernel=rr(size(estimate))<5;
    mypos=floor(size(estimate)/2);
    tmp(mypos(1),mypos(2))=1;
    mymask=convolve(tmp,mykernel) > 0.5; %just one point in the middle
%     myillu_mask=repmat(mymask,[1 1 NumIm]);
for i=1:NumIm
    myillu_mask{i}=mymask;
%     myillu_mask{i}=fft2rft(mymask);
end
end

%%   To not rerun the random generators

%disableCuda();
h=obj*0;h(2,2)=1;h=gaussf(h);
img=sqrt(prod(size(obj)))*real(ift(ft(obj) .* ft(h)));
%oimg=convolve(obj,h);

global myim;myim={};myim{1}=img; myim{2}=1.3*img+1;
global otfrep;otfrep={};otfrep{1}=rft(h);
global lambdaPenalty;lambdaPenalty=1.0;
global DeconvMethod;DeconvMethod='LeastSqr';
%global DeconvMethod;DeconvMethod='Poisson';
global NegPenalty;NegPenalty='NONE';
%global NegPenalty;NegPenalty='NegSqr';
%global RegularisationMethod;RegularisationMethod='GR';
%global RegularisationMethod;RegularisationMethod='GS';
%global RegularisationMethod;RegularisationMethod='AR';
%global RegularisationMethod;RegularisationMethod='TV';
global RegularisationMethod;RegularisationMethod='NONE';
global DeconvVariance;
DeconvVariances=[];
global BetaVals;
BetaVals=[1 1 1];
%global DeconvMask;DeconvMask=xx(sX,sY)>0;  % only data in this mask will be evaluated
global DeconvMask;DeconvMask=[];  % only data in this mask will be evaluated
%global myillu;myillu{1}=illu;  % only data in this mask will be evaluated
global myillu;myillu={}; % myillu{1}=illu(:,:,0);myillu{2}=illu(:,:,1); 
global NormFac;NormFac=1.0;   % Normalisation factor
global aRecon;aRecon=objestimate;


% global PupilInterpolators;
% PupilInterpolators={};
global RegularisationParameters;
ToReg=1;  % 0 is object, 1 means illu, 2 means otf
% RegularisationParameters=ParseRegularisation({{'ProjPupil',[488,1.4,1.518],[100 100 100]}},ToReg); % will set RegularisationParameters(14,1)=1 % case 'ProjPupil' where only the 2D pupil is estimated
RegularisationParameters=ParseRegularisation(MyReg,ToReg);
% RegularisationParameters=ParseRegularisation({'ForcePos',[]},ToReg);

global ToEstimate;ToEstimate=1;   % 0 is object (with or without known illu), 1 is illu, 2 is OTF
AssignFunctions(RegularisationParameters,2)  % here: 0 is object, 1 is object with known illu, 2 is illum, 3 is OTF
%myVec=double(reshape(estimate,prod(size(estimate)))); % object estimate
% global OTFmask;
% OTFmask=PupilInterpolators.Mask;

% mysize=size(PupilInterpolators.indexList2D,2);
mysize=[10 10];
% VecOTF=(rand(2*mysize,1)-0.5)*2;  % for real and imaginary part
% myVec=transpose(VecOTF);
% myVec=double(reshape(estimate(OTFmask),prod(size(estimate)))); % object estimate

% Old version:
myVec=double(reshape(estimate,prod(size(estimate))));

global myillu_sumcond;
global my_sumcond
if (~ useSumCond)
    myillu_sumcond={};
    my_sumcond={};
    myVec=[myVec myVec];
else
    myillu_sumcond={NumIm};
    my_sumcond={NumIm};
end

% myillu_sumcond={};


%%
myVec=2.0+myVec*0;
%%
[err,grad]=GenericErrorAndDeriv(myVec);  % to get the analytical gradient solution
%%   Do all the cuda conversions here to test the cuda variant of the algorithm
if (useCuda)
initCuda();
myVec=cuda(myVec);
grad=cuda(grad);
for n=1:numel(myillu)
    myillu{n}=cuda(myillu{n});
    myim{n}=cuda(myim{n});
end
otfrep{1}=rft(fftshift(cuda(h)));
DeconvMask=ConditionalCudaConvert(DeconvMask,useCuda);
aRecon=ConditionalCudaConvert(aRecon,useCuda);

[err,grad]=GenericErrorAndDeriv(myVec);
end
%%
clear mygrad;
eps = 1e-3;
fprintf('Testing Gradient direction total: %g\n',size(myVec,2));
for d=1:size(myVec,2)
    fprintf('%d ',d);
    UnitD = myVec*0;
    UnitD(d) = 1;
    mygrad(d) = (GenericErrorAndDeriv(myVec+(eps * UnitD)) - err) / eps;
    if mod(d,40)==0
        fprintf('\n');
    end
end
fprintf('\n');

%grad= reshape(dip_image(grad','single'),size(img));
%mygrad=reshape(dip_image(mygrad','single'),size(img));
global savedATF;
global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF

% AssignToGlobal(ConvertInputToModel(grad)); % changes otfrep{1} and SavedATF
% grad=savedATF;
if(1)
    if (useSumCond)
        grad= reshape(dip_image(grad','single'),size(img)); %bug if rft and mask
    else
        grad= reshape(dip_image(grad','single'),[size(img) 2]); %bug if rft and mask
    end
else
    grad2=extract(dip_image(grad','single'),[size(img,1)]); %because of the mask, length is too short
    grad2=reshape(dip_image(grad2','single'),[ceil(sqrt(prod(size(grad2)))) ceil(sqrt(prod(size(grad2))))]);
    % grad2= reshape(dip_image(grad','single'),[ceil(sqrt(prod(size(grad)))) ceil(sqrt(prod(size(grad))))]); %didnt work
    grad3=rft2fft(grad2); %should now be the same size as mygrad
end

% AssignToGlobal(ConvertInputToModel(mygrad)); % changes otfrep{1} and SavedATF
% mygrad=savedATF;
    if (useSumCond)
        mygrad=reshape(dip_image(mygrad','single'),size(img));
    else
        mygrad=reshape(dip_image(mygrad','single'),[size(img) 2]);
    end

if isa(grad,'cuda')
    mygrad=cuda(mygrad);
end

todisplay=cat(3,grad,mygrad)

relerror = (mygrad - grad') ./ max(abs(grad'));
fprintf('Max Error :%g\n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges
if (useSumCond)
    fprintf('Max Center Error :%g\n',max(abs(relerror(1:end-1,1:end-1))))   % Problems are caused by the hessian operator on finite arrays at the edges
else
    fprintf('Max Center Error :%g\n',max(abs(relerror(1:end-1,1:end-1,:))))   % Problems are caused by the hessian operator on finite arrays at the edges
end