%03.03.2014: Test new deconv toolbox OK. Aurelie
if exist('cuda_cuda')
    disableCuda();
end

sX=10;
sY=10;
NumIm=2;
useCuda=0;
MyObjReg={'ForcePhase',[];'NormMeasSumSqr',[]};  % ;'TV',[0.1,0]; 'ER',[1 1]; 

rng(1); % initialize the random generator with the same seed always
r1=rand(sX,sY);
rng(2); % initialize the random generator with the same seed always
% obj = dip_image(r1);
r2=rand(sX,sY);
obj=dip_image(r1+1i*r2);

%%estimate=mean(img)+img*0;  % Careful: This yields rubbish with estimating Regularisations
%%estimate=obj*1.1;
%%estimate=newim(9,9);estimate(4,4)=0.5;

%obj=newim(9,9);obj(4,4)=1;
rng(3); % initialize the raqdom generator with the same seed always
r1=rand(sX,sY);
%rng(4); % initialize the raqdom generator with the same seed always
%r2=rand(sX,sY);
estimate=dip_image(r1);
rng(5); % initialize the raqdom generator with the same seed always
r1=rand(sX,sY);
rng(6); % initialize the random generator with the same seed always
% objestimate = dip_image(r1);
r2=rand(sX,sY);
objestimate=dip_image(r1+1i*r2);
rng(7); % initialize the random generator with the same seed always
illu=dip_image(rand(sX,sY,NumIm));   % an illumination distribution

%%   To not rerun the random generators

atf = newim(size(obj)); 
mymask1 = (abs(xx(atf))<=1); mymask2 = (abs(yy(atf))<=1);
mymask = mymask1.*mymask2;
atf(mymask) = 1;
img=sqrt(prod(size(obj)))*(ift(ft(obj) .* atf));
% h=obj*0;h(2,2)=1;h=gaussf(h);
% img=sqrt(prod(size(obj)))*real(ift(ft(obj) .* ft(h)));
%oimg=convolve(obj,h);

global myim;myim={};myim{1}=img;% myim{2}=img;

global otfrep;otfrep={};otfrep{1}=atf; % rft(abssqr(atf));
global lambdaPenalty;lambdaPenalty=1.0;
global DeconvMethod;DeconvMethod='LeastSqr';
%global DeconvMethod;DeconvMethod='Poisson';
global NegPenalty;NegPenalty='NONE';
global DeconvVariance; DeconvVariances=[];
global BetaVals;
BetaVals=[1 1 1];
global DeconvMask;DeconvMask=[]; %xx(sX,sY)>0;  % only data in this mask will be evaluated
%global DeconvMask;DeconvMask=[];  % only data in this mask will be evaluated
%global myillu;myillu{1}=illu;  % only data in this mask will be evaluated
global myillu;myillu={};myillu{1}=illu(:,:,0);myillu{2}=illu(:,:,1);  % only data in this mask will be evaluated
global NormFac;NormFac=1.0;   % Normalisation factor
global ToEstimate;ToEstimate=0;   % 0 is object, 1 is illu
ToReg=0;  % 0: Object, 1 means illu
global aRecon;aRecon=objestimate;   % 0 is object, 1 is illu
global ComplexPSF; ComplexPSF=1; %check that this is correct. Aurelie
global RegularisationParameters;
RegularisationParameters=ParseRegularisation(MyObjReg,ToReg);
AssignFunctions(RegularisationParameters,1)
myVec=double(reshape(estimate,prod(size(estimate))));
global measSums;
global measSumsSqr;
global measSumSqr;
v=1;
measSums{v}=sum(myim{v},DeconvMask);
measSumsSqr{v}=sqrt(sum(real(myim{v} .* conj(myim{v})),DeconvMask));
measSumSqr=measSumsSqr{1};
% myVec = [real(myVec), imag(myVec)];
%%
[err,grad]=GenericErrorAndDeriv(myVec);
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
    DeconvMask=cuda(DeconvMask);
    aRecon=cuda(aRecon);
    
    [err,grad]=GenericErrorAndDeriv(myVec);
end

%%
clear mygrad;
eps = 2e-4;
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

% grad(1:end/2) = grad(1:end/2)+1i*grad(end/2+1:end);
% grad(end/2+1:end) = [];
% mygrad(1:end/2) = mygrad(1:end/2)+1i*mygrad(end/2+1:end);
% mygrad(end/2+1:end) = [];
grad= reshape(dip_image(grad','scomplex'),size(img));
mygrad=reshape(dip_image(mygrad','scomplex'),size(img));

if isa(grad,'cuda')
    mygrad=cuda(mygrad);
end

cat(3,grad,mygrad)

relerror = (mygrad - grad') ./ mean(abs(grad));
fprintf('Max Error :%g\n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges
fprintf('Max Center Error :%g\n',max(abs(relerror(1:end-1,1:end-1))))   % Problems are caused by the hessian operator on finite arrays at the edges
