cd y:\MATLAB\Toolboxes\Deconv\DeconvolutionTests\
useCuda=0;
a=load('D:\Data\BlindDeconvData_Autoquant\s.mat');
img=squeeze(a.s);

%%
if useCuda
    disableCuda();
end
scaleX=64.5;
scaleY=scaleX;
scaleZ=200;
NA=1.4;
lambda=520;
ri=1.52;

%% Intensity-based estimation

% [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[5 150 0 0 0],'Poisson','RLL',{{},{},{}},[1,1,1],[0 0 0],[],useCuda); 
% Intensity-based estimation
h=aberratedPSF(0,scaleX,lambda,NA,{'ri',ri;'lambdaEm',lambda;'na',NA;'scaleZ',scaleZ;'sX',size(img,1);'sY',size(img,2);'sZ',size(img,3);});
[myDeconvResInt,dum,resPSFInt]=GenericDeconvolution(img,h,[15 150 50 0 50],'LeastSqr','Blbfgs',{{'ForcePos',1;},{},{'NegSqr',0.05;}},[1,1,1],[0 0 0],[],useCuda);

%% Amplitude-based blind PSF estimation:

smallimg=extract(img,[128 128 16],[146 146 22]);
% [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,[],[16 5 5 0 40],'LeastSqr','Blbfgs',{{'ForcePos',1},{},{'ProjPupil',[lambda,NA,ri],[scaleX scaleY scaleZ]}},[scaleX scaleY scaleZ],[0 0 9],[],useCuda);
NA_deconv=1.4;
scaleZ=85;  % There is something wrong with the Z scaling in this data
[myDeconvRes,dum,resPSF]=GenericDeconvolution(smallimg,[],[16 5 5 0 10],'LeastSqr','Blbfgs',{{'ForcePos',1},{},{'ProjPupil',[lambda,NA_deconv,ri],[scaleX scaleY scaleZ]}},[scaleX scaleY scaleZ],[0 0 9],[],useCuda);

res3=myDeconvRes;
psf3=resPSF;
global savedASF
atf3=savedASF;
%myDeconvRes2=GenericDeconvolution(img,resPSF,150,'LeastSqr','RLL',{{},{},{}},[scaleX scaleY scaleZ],[0 0 9],[],useCuda);

%% Try a supersampled deconv
 resampledPSF=real(ift(extract(ft(resPSF{1}),[size(resPSF{1}).*[2 2 1]])));
[myDeconvResInt,dum,resPSFInt]=GenericDeconvolution(img,resampledPSF,[15 150 0 0 0],'LeastSqr','lbfgs',{{'ForcePos',1;'Resample',[2 2 1]},{},{}},[1,1,1],[0 0 0],[],useCuda);

%%
[myDeconvR1,dum,resPSF1]=GenericDeconvolution(img,resPSF,150,'LeastSqr','lbfgs',{{'ForcePos',1;'Resample',[0.5 0.5 1]},{},{}},[1,1,1],[0 0 0],[],useCuda);

