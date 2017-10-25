%03.03.2014: Test new deconv toolbox OK. Aurelie
% 02.06.2015 : Introduced complex numbers for the object estimate to test the performance. Seems to work, but there were errors in the comparison in this file for complex numbers

% This function simulates a small image and test numerically whether the gradiests as estimated by GenericErrorAndDrivative correspond to what is expected from the error value return of that function
%clear all
%disableCuda();
useCuda=0; %disableCuda();
sX=10; %10 and 50: checked
sY=10;
sZ=4;
NumIm=2;

obj=dip_image(rand(sX,sY,sZ));

%%estimate=mean(img)+img*0;  % Careful: This yields rubbish with estimating Regularisations
%%estimate=obj*1.1;
%%estimate=newim(9,9);estimate(4,4)=0.5;

%obj=newim(9,9);obj(4,4)=1;
% estimate=dip_image(rand(sX,sY)+i*rand(sX,sY));  % new OTF estimate (everywhere, not only in the mask)
objestimate=dip_image(rand(sX,sY,sZ));  % current (old) reconstruction estimate
illu=dip_image(rand(sX,sY,sZ,NumIm));   % an illumination distribution
%estimate=dip_image(rand(sX,sY,sZ));
estimate=dip_image(rand(sX,sY,sZ)+i*rand(sX,sY,sZ));


%%   To not rerun the random generators

%disableCuda();
h=obj*0;h(4,4,1)=1;h=gaussf(h);
img=sqrt(prod(size(obj)))*real(ift(ft(obj) .* ft(h)));
%oimg=convolve(obj,h);

global myim;myim={};myim{1}=img; myim{2}=img;
global otfrep;otfrep={};otfrep{1}=rft(h);
global lambdaPenalty;lambdaPenalty=1.0;
global DeconvMethod;DeconvMethod='LeastSqr';
%global DeconvMethod;DeconvMethod='Poisson';
global NegPenalty;NegPenalty='NONE';
%global NegPenalty;NegPenalty='NegSqr';
%global RegularisationMethod;RegularisationMethod='GR';
global RegularisationMethod;RegularisationMethod='ER';
%global RegularisationMethod;RegularisationMethod='AR';
%global RegularisationMethod;RegularisationMethod='TV';
% global RegularisationMethod;RegularisationMethod='NONE';
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
ToReg=0;  % 0 is object, 1 means illu, 2 means otf
% RegularisationParameters=ParseRegularisation({{'ProjPupil',[488,1.4,1.518],[100 100 100]}},ToReg); % will set RegularisationParameters(14,1)=1 % case 'ProjPupil' where only the 2D pupil is estimated
Regu={{'ForcePos',[];'GS',1e-4},{'ForcePos',[]}};
Regu={{'GS',1e-4},{}};
RegularisationParameters=ParseRegularisation(Regu,ToReg);
global ToEstimate;ToEstimate=1;   % 0 is object (with or without known illu), 1 is illu, 2 is OTF
global ComplexPSF; ComplexPSF=0; %check that this is correct. Aurelie
AssignFunctions(RegularisationParameters,0)  % here: 0 is object, 1 is object with known illu, 2 is illum, 3 is OTF
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
myillu_sumcond={NumIm};
% myillu_sumcond={};

global my_sumcond
my_sumcond={NumIm};

%%
% [err,grad]=GenericErrorAndDeriv(myVec);  % to get the analytical gradient solution
[err,grad]=Regularize(myVec,BetaVals); %input should be estimate?

%%   Do all the cuda conversions here to test the cuda variant of the algorithm
if (1)
initCuda();
myVec=cuda(myVec);
[err,grad]=Regularize(myVec,BetaVals); %input should be estimate?

% grad=cuda(grad);
% for n=1:numel(myillu)
%     myillu{n}=cuda(myillu{n});
%     myim{n}=cuda(myim{n});
% end
% otfrep{1}=rft(fftshift(cuda(h)));
% % DeconvMask=cuda(DeconvMask);
% aRecon=cuda(aRecon);
% 
% [err,grad]=GenericErrorAndDeriv(myVec);
end
%%
clear mygrad;
eps = 1e-3;
% myVec=estimate; %To avoid renaming everywhere
for d=1:size(myVec,2)
    UnitD = myVec*0;
    UnitD(d) = 1;
%     mygrad(d) = (GenericErrorAndDeriv(myVec+(eps * UnitD)) - err) / eps;
    mygrad(d) = (Regularize(myVec+(eps * UnitD),BetaVals) - err) / eps;
    mygrad(d) = mygrad(d)+ i*(Regularize(myVec+(eps * i*UnitD),BetaVals) - err) / eps;
end

%grad= reshape(dip_image(grad','single'),size(img));
%mygrad=reshape(dip_image(mygrad','single'),size(img));
% global savedATF;
% global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
% global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF

% AssignToGlobal(ConvertInputToModel(grad)); % changes otfrep{1} and SavedATF
% grad=savedATF;

% grad= reshape(dip_image(grad','single'),size(img)); %bug if rft and mask
grad= reshape(dip_image(transpose(grad),'scomplex'),size(img)); %bug if rft and mask

% AssignToGlobal(ConvertInputToModel(mygrad)); % changes otfrep{1} and SavedATF
% mygrad=savedATF;
%mygrad=reshape(dip_image(mygrad','single'),size(img));
mygrad=reshape(dip_image(transpose(mygrad),'scomplex'),size(img));

if isa(grad,'cuda')
    mygrad=cuda(mygrad);
end

cat(4,grad,mygrad)

relerror = (mygrad - grad) ./ max(abs(grad));
fprintf('Max Error :%g\n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges
fprintf('Max Center Error :%g\n',max(abs(relerror(1:end-1,1:end-1,:))))   % Problems are caused by the hessian operator on finite arrays at the edges
