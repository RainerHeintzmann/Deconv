%clear all
disableCuda();

sX=10;
sY=10;
NumIm=2;

obj=dip_image(rand(sX,sY));

%%estimate=mean(img)+img*0;  % Careful: This yields rubbish with estimating Regularisations
%%estimate=obj*1.1;
%%estimate=newim(9,9);estimate(4,4)=0.5;

%obj=newim(9,9);obj(4,4)=1;
estimate=dip_image(rand(sX,sY));
objestimate=dip_image(rand(sX,sY));
illu=dip_image(rand(sX,sY,NumIm));   % an illumination distribution

%%   To not rerun the random generators

disableCuda();
h=obj*0;h(2,2)=1;h=gaussf(h);
img=sqrt(prod(size(obj)))*real(ift(ft(obj) .* ft(h)));
%oimg=convolve(obj,h);

global myim;myim={};myim{1}=img;myim{2}=img;
global otfrep;otfrep={};otfrep{1}=ft(h);
global lambdaPenalty;lambdaPenalty=10.0;
global DeconvMethod;DeconvMethod='LeastSqr';
%global DeconvMethod;DeconvMethod='Poisson';
global NegPenalty;NegPenalty='NONE';
%global NegPenalty;NegPenalty='NegSqr';
%global RegularisationMethod;RegularisationMethod='GR';
%global RegularisationMethod;RegularisationMethod='AR';
%global RegularisationMethod;RegularisationMethod='TV';
global RegularisationMethod;RegularisationMethod='NONE';
global DeconvVariance;
DeconvVariances=[];
global BetaVals;
BetaVals=[1 1 1];
global DeconvMask;DeconvMask=xx(sX,sY)>0;  % only data in this mask will be evaluated
%global DeconvMask;DeconvMask=[];  % only data in this mask will be evaluated
%global myillu;myillu{1}=illu;  % only data in this mask will be evaluated
global myillu;myillu={};myillu{1}=illu(:,:,0);myillu{2}=illu(:,:,1);  % only data in this mask will be evaluated
global NormFac;NormFac=1.0;   % Normalisation factor
global ToEstimate;ToEstimate=1;   % 0 is object, 1 is illu
global aRecon;aRecon=objestimate;   % 0 is object, 1 is illu

myVec=double(reshape(estimate,prod(size(estimate))));
[err,grad]=GenericErrorAndDeriv(myVec);
%%   Do all the cuda conversions here to test the cuda variant of the algorithm
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

%%
clear mygrad;
eps = 1e-3;
for d=1:size(myVec,2)
    UnitD = myVec*0;
    UnitD(d) = 1;
    mygrad(d) = (GenericErrorAndDeriv(myVec+(eps * UnitD)) - err) / eps;
end

grad= reshape(dip_image(grad','single'),size(img));
mygrad=reshape(dip_image(mygrad','single'),size(img));

if isa(grad,'cuda')
    mygrad=cuda(mygrad);
end

cat(3,grad,mygrad)

relerror = (mygrad - grad') ./ max(abs(grad'));
fprintf('Max Error :%g\n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges
fprintf('Max Center Error :%g\n',max(abs(relerror(1:end-1,1:end-1))))   % Problems are caused by the hessian operator on finite arrays at the edges
