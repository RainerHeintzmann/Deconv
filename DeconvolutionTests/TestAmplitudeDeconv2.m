%%  Now try this in 3D!
addpath('y:\MATLAB\Toolboxes\holography\');
addpath('y:\MATLAB\Toolboxes\Objects\');
disableCuda();
% enableCuda();
mysize3d=[100 90 100];
imatrix=IterateCoefficients(40,6,mysize3d(3),10,500);  % 40 subpixel subdivisions, 2*10+1 kernelsize, 20 pixel bordersize in all directions, 1500 iterations
%
aproj=extract(readim,mysize3d(1:2));

%% make a 3D PSF by filling a part of a spherical shell in Fourier space
%k0=40; kxymax=22;
maxphotons=1000;
lambda=515; pixelsize=[100 100 100]; NA=0.8;refractiveIndex=1.0;
tic
%[indexList2D,fullIndex3D,factorList]=FillProjSpherePrepare(size3d,k0,kxymax,imatrix);  % old syntax
[indexList2D,fullIndex3D,factorList,aMask]=FillProjSpherePrepare(mysize3d,lambda,pixelsize,NA,imatrix,refractiveIndex);
toc
famp=newim(mysize3d,'scomplex'); % empty 3d complex array for Fourier-space
mypupil=aproj*0+1;
%mypupil=readim+i*readim('orka');
%mypupil=extract(mypupil,size3d(1:2));
atf=FillProjSphere(famp,mypupil,indexList2D,fullIndex3D,factorList);  % Just a circular wave or a complicated phase distribution
pupil=ProjSphere2D(repmat(atf,[1 1 1 1]),indexList2D,fullIndex3D,factorList)  % estimates the 2D pupil from a given 3D distribution
toc
asf=ift(atf);
toc
%%
if (0)
    obj=Gen3dObj(size(asf),100,100,100,2,2);
else
    obj=newim(size(asf));
end
numPts=50; %0;
obj(floor(prod(size(obj))*rand(numPts,1))-1)=30;

%theillu=exp(-i*2*pi*zz(obj)/(lambda/pixelsize(3)));
theillu=exp(-i*2*pi*zz(obj)/(lambda/pixelsize(3)));

myamp=ift(atf .*ft(obj .* theillu));
myint=abssqr(myamp);
%%   Do reconstructions with the 3D backpropagated amplitude distribution as an input.
usecuda=1;
% First we ignore that the object was illuminated, as this can also be interpreted as part of the object phase.
global myillu
myillu=[];
resamp=GenericDeconvolution(myamp,asf,100,'LeastSqr','lbfgs',{'Complex',[]},[1 1 1],[0 0 0],[],usecuda)  % Looks OK, but is still only very high frequency
% result is not so bad.  Why is the ring shifted to even higher frequencies??

resampC=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'Complex',[];'GS',100;'StartImg',abs(resamp);'Illumination',{theillu}},[1 1 1],[0 0 0],[],usecuda)
resampC=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'Complex',[];'GS',10;'StartImg',resampC;'Illumination',{theillu}},[1 1 1],[0 0 0],[],usecuda)
resampC=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'Complex',[];'CO',-10;'StartImg',abs(resamp);'Illumination',{theillu}},[1 1 1],[0 0 0],[],usecuda)
cat(4,myint,gaussf(obj,[1 1 2]),abssqr(resamp),abssqr(resampC))
%resamp=GenericDeconvolution(myamp,asf,100,'LeastSqr','lbfgs',{'Complex',[];'GS',100;'Illumination',{theillu}},[1 1 1],[0 0 0],[],usecuda) % The regularisation does not really help.
% How about a phase-only regularizer?

%resamp2=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'TV',[1,0];'ForcePos',[]},[1 1 10],[0 0 0],[],usecuda)
%resamp2=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'ForcePos',[];'StartImg',abssqr(resamp);},[1 1 1],[0 0 0],[],usecuda)
myillu=[];
resamp2=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'ForcePos',[];'StartImg',abssqr(resamp);'GS',1},[1 1 1],[0 0 0],[],usecuda)
resamp3=GenericDeconvolution(myamp,asf,20,'LeastSqr','lbfgs',{'ForcePos',[];'StartImg',abssqr(resamp);'GS',10;'Illumination',{theillu}},[1 1 1],[0 0 0],[],usecuda)

resamp4=GenericDeconvolution(myamp,asf,100,'LeastSqr','lbfgs',{'ForcePos',[];'GS',100;'Illumination',{theillu}},[1 1 1],[0 0 0],[],usecuda)
%recon=GenericDeconvolution(myint,asf,20,'LeastSqr','lbfgs',{'IntensityData',[];},[1 1 1],[0 0 0],[],usecuda)

% Now try with shifted ATF and no illumination.  -> even worse!
resamp5=GenericDeconvolution(myamp*conj(theillu),asf,20,'LeastSqr','lbfgs',{'ForcePos',[];'StartImg',abssqr(resamp);'GS',100},[1 1 1],[0 0 0],[],usecuda)
resamp6=GenericDeconvolution(myamp*conj(theillu),asf,20,'LeastSqr','lbfgs',{'ForcePos',[];'GS',1},[1 1 1],[0 0 0],[],usecuda)

cat(4,myint,gaussf(obj,[1 1 2]),abssqr(resamp),abssqr(resamp2),abssqr(resamp3),gaussf(obj,[1 1 2]))
