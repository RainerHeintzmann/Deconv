% This function simulates a small image and test numerically whether the gradiests as estimated by GenericErrorAndDrivative correspond to what is expected from the error value return of that function
%clear all
%disableCuda();

sX=10;
sY=10; 
NumIm=2;
rng(0);
obj=dip_image(rand(sX,sY));

%%estimate=mean(img)+img*0;  % Careful: This yields rubbish with estimating Regularisations
%%estimate=obj*1.1;
%%estimate=newim(9,9);estimate(4,4)=0.5;

%obj=newim(9,9);obj(4,4)=1;
% estimate=dip_image(rand(sX,sY)+i*rand(sX,sY));  % new OTF estimate (everywhere, not only in the mask)
rng(1);
objestimate=dip_image(rand(sX,sY));  % current (old) reconstruction estimate
rng(2);
illu=dip_image(rand(sX,sY,NumIm));   % an illumination distribution

%%   To not rerun the random generators

%disableCuda();
h=obj*0;h(3,2)=1;h=gaussf(h);
img=sqrt(prod(size(obj)))*real(ift(ft(obj) .* ft(h)));
%oimg=convolve(obj,h);

global myim;myim={};myim{1}=img; % myim{2}=img;
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


global PupilInterpolators;
PupilInterpolators={};
global RegularisationParameters;
ToReg=2;  % 0 is object, 1 means illu, 2 means otf
RegularisationParameters=ParseRegularisation({{'ProjPupil',[488,1.4,1.518],[100 100 100]}},ToReg); % will set RegularisationParameters(14,1)=1 % case 'ProjPupil' where only the 2D pupil is estimated
global ToEstimate;ToEstimate=2;   % 0 is object (with or without known illu), 1 is illu, 2 is OTF
AssignFunctions(RegularisationParameters,3)  % here: 0 is object, 1 is object with known illu, 2 is illum, 3 is OTF
%myVec=double(reshape(estimate,prod(size(estimate)))); % object estimate
global OTFmask;
OTFmask=PupilInterpolators.Mask;

mysize=size(PupilInterpolators.indexList2D,2);
rng(3);
VecOTF=(rand(2*mysize,1)-0.5)*2;  % for real and imaginary part
myVec=transpose(VecOTF);
myVec(1:mysize)=1;myVec(mysize+1:end)=0;
% myVec=double(reshape(estimate(OTFmask),prod(size(estimate)))); % object estimate

%%
[err,grad]=GenericErrorAndDeriv(myVec);  % to get the analytical gradient solution
%%   Do all the cuda conversions here to test the cuda variant of the algorithm
if (0)
initCuda();
myVec=cuda(myVec);
grad=cuda(grad);
for n=1:numel(myillu)
    myillu{n}=cuda(myillu{n});
    myim{n}=cuda(myim{n});
end
otfrep{1}=rft(fftshift(cuda(h)));
DeconvMask=cuda(DeconvMask);
aRecon=cuda(aRecon);

[err,grad]=GenericErrorAndDeriv(myVec);
end
%%
clear mygrad;
eps = 2e-3;
for d=1:size(myVec,2)
    UnitD = myVec*0;
    UnitD(d) = 1;
    mygrad(d) = (GenericErrorAndDeriv(myVec+(eps * UnitD)) - err) / eps;
end

%grad= reshape(dip_image(grad','single'),size(img));
%mygrad=reshape(dip_image(mygrad','single'),size(img));
global savedATF;
global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF

AssignToGlobal(ConvertInputToModel(grad)); % changes otfrep{1} and SavedATF
grad=savedATF;
%grad= reshape(dip_image(grad','single'),size(img));
AssignToGlobal(ConvertInputToModel(mygrad)); % changes otfrep{1} and SavedATF
mygrad=savedATF;
%mygrad=reshape(dip_image(mygrad','single'),size(img));

if isa(grad,'cuda')
    mygrad=cuda(mygrad);
end

cat(3,grad,mygrad)

relerror = (mygrad - grad) ./ mean(abs(grad));
fprintf('Max Error :%g\n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges
fprintf('Max Center Error :%g\n',max(abs(relerror(1:end-1,1:end-1))))   % Problems are caused by the hessian operator on finite arrays at the edges
