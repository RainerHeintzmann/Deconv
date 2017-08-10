% The aim is to leave one out and evaluate the error of the left out pixels in dependence of the regularisation factor lambda.
% This allows to estimate lambda
obj=readim('Y:\MATLAB\images\resolution_coarse')+10;
nphoton=100;
ImageParam=struct('Sampling',[30 30 200],'Size',[size(obj,1) size(obj,2) 1]);
PSFParam=struct('NA',1.4,'n',1.518,'MinOtf',1.2e-3,'lambdaEm',520);

[h,amp]=GenericPSFSim(ImageParam,PSFParam,'VolumeShell');  % has problems with undersampling in Z
otf=ft(squeeze(h));
pimg=real(ift(ft(obj) .* otf));
pimg=pimg./max(pimg)*nphoton;
%nimg=noise(pimg,'Gaussian',100)
nimg=noise(pimg,'Poisson')
%

%myDeconv=GenericDeconvolution(nimg,h,150,'LeastSqr','lbfgs',{'ForcePos',1},[1 1 1],[0 0 0],[],useCuda); 
myDeconv=GenericDeconvolution(nimg,h,150,'Poisson','lbfgs',{'ForcePos',1},[1 1 1],[0 0 0],[],useCuda); 
mydiff=rift(Recons)-nimg;
WrongPixels=sum(mydiff(~myMask).^2);
st=cat(1,img{1},myDeconv)

myMask=nimg*0+1;
myMask(floor(rand(1,floor(prod(size(myMask))*0.1))*prod(size(myMask))))=0;
myMask=myMask>0.5;

mylambda=0.02;
%myDeconvR=GenericDeconvolution(nimg,h,50,'LeastSqr','lbfgs',{'TV',[mylambda,1e-6];'DeconvMask',myMask},[1 1 1],[0 0 0],[],useCuda); % 'ForcePos',1;
myDeconvR=GenericDeconvolution(nimg,h,100,'Poisson','lbfgs',{'ForcePos',1;'GR',mylambda;'DeconvMask',myMask},[1 1 1],[0 0 0],[],useCuda); % 'ForcePos',1;
global FinalErrNorm;
global ConvertModelToVec;
global NormFac;
global DeconvMask;
EstimationErr=[];
myInd=1;
for mylambda=0:1e-6:1e-5
    myDeconvR=GenericDeconvolution(nimg,h,50,'Poisson','lbfgs',{'NormFac',1.0;'ForcePos',1;'GS',mylambda;'DeconvMask',myMask;'Reuse',1},[1 1 1],[0 0 0],[],useCuda); % 'ForcePos',1;
    dipshow(1,myDeconvR);drawnow();
    DeconvMask=~myMask;
    myRes=ConvertModelToVec(myDeconvR);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation
    [myErr myDeriv]=GenericErrorAndDeriv(myRes);
    EstimationErr(myInd)=myErr / NormFac
    myInd=myInd+1;
end
mydiffR=Recons-nimg;
WrongPixelsR=sum(mydiffR(~myMask).^2);

