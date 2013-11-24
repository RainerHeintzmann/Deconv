% First test with just a simple amplitude psf in 2D
% obj=readim+i*readim('orka');
%disableCuda();
maxphotons=1000;
obj=readim;
atf=(rr<30) * exp(i*2*pi*(100*rr('freq').^2)); % Fresnel defocus
asf=ift(atf);
dat=ift(atf .* ft(obj));

myint=abssqr(dat);

dat=dat/sqrt(max(myint)*maxphotons);
myint=noise(myint/max(myint)*maxphotons,'poisson');


usecuda=0;
% deconvolutions using measured amplitudes
res=GenericDeconvolution(dat,asf,100,'LeastSqr','lbfgs',{'Complex',[]},[1 1],[0 0],[],usecuda)
resTV=GenericDeconvolution(dat,asf,100,'LeastSqr','lbfgs',{'TV',[0.001,0];'Complex',[]},[1 1],[0 0],[],usecuda)
%%
% try to get arbitrary amplitudes from intensity data
resA=GenericDeconvolution(myint,asf,100,'LeastSqr','lbfgs',{'Complex',[];'IntensityData',[];},[1 1],[0 0],[],usecuda);
resAbsAmp=abs(resA)
%resIA=GenericDeconvolution(myint,asf,10,'LeastSqr','lbfgs',{'TV',[0.1,100.0];'Complex',[];'IntensityData',[];},[1 1],[0 0],[],0)
% amplitudes from intensity data but with real scattering constraints
resAReal=GenericDeconvolution(myint,asf,100,'LeastSqr','lbfgs',{'IntensityData',[];},[1 1],[0 0],[],usecuda)

resInt=GenericDeconvolution(myint,abssqr(asf),100,'LeastSqr','lbfgs',[],[1 1],[0 0],[],usecuda); % standard intensity deconv

%%  Now try this in 3D!
addpath('y:\MATLAB\Toolboxes\holography\');
addpath('y:\MATLAB\Toolboxes\Objects\');
%disableCuda();
% enableCuda();
size3d=[100 100 100];
imatrix=IterateCoefficients(40,6,size3d(3),10,500);  % 40 subpixel subdivisions, 2*10+1 kernelsize, 20 pixel bordersize in all directions, 1500 iterations
%
aproj=readim;

%% make a 3D PSF by filling a part of a spherical shell in Fourier space
k0=40; kxymax=22;
famp=newim(size3d,'scomplex'); % empty 3d complex array for Fourier-space
tic
[indexList2D,fullIndex3D,factorList]=FillProjSpherePrepare(size3d,k0,kxymax,imatrix);
toc
atf=FillProjSphere(famp,aproj*0+1,indexList2D,fullIndex3D,factorList);  % Just a circular wave
toc
asf=ift(atf);
toc
%%
obj=Gen3dObj(size(asf),100,100,100,2,2);
numPts=500;
obj(floor(prod(size(obj))*rand(numPts,1))-1)=30;

myamp=ift(atf .*ft(obj));
myint=abssqr(myamp);
%%
usecuda=0;
resamp=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'Complex',[]},[1 1 1],[0 0 0],[],usecuda)

%resamp2=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'TV',[1,0];'ForcePos',[]},[1 1 10],[0 0 0],[],usecuda)
%resamp2=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'ForcePos',[];'StartImg',abssqr(resamp);},[1 1 1],[0 0 0],[],usecuda)
resamp2=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'ForcePos',[];'StartImg',abssqr(resamp);'GS',1},[1 1 1],[0 0 0],[],usecuda)
 
recon=GenericDeconvolution(myint,asf,20,'LeastSqr','lbfgs',{'IntensityData',[];},[1 1 1],[0 0 0],[],usecuda)

cat(4,myint,gaussf(obj,[1 1 2]),abssqr(resamp),abssqr(resamp2))
