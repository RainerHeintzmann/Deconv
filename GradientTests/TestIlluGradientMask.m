% This function simulates a small image and test numerically whether the gradiests as estimated by GenericErrorAndDrivative correspond to what is expected from the error value return of that function
%clear all
%disableCuda();
useCuda=0; % disableCuda();
sX=102; %10 and 50: checked
sY=100;
NumIm=4;
if (0)
    MyIlluReg={};  % ;'TV',[0.1,0]
else
    MyIlluReg={'ForcePos',[]};  % Why does this not work?
end

obj=dip_image(rand(sX,sY));

%%estimate=mean(img)+img*0;  % Careful: This yields rubbish with estimating Regularisations
%%estimate=obj*1.1;
%%estimate=newim(9,9);estimate(4,4)=0.5;

%obj=newim(9,9);obj(4,4)=1;
% estimate=dip_image(rand(sX,sY)+i*rand(sX,sY));  % new OTF estimate (everywhere, not only in the mask)
objestimate=10*dip_image(rand(sX,sY));  % current (old) reconstruction estimate
illu=dip_image(rand(sX,sY,NumIm));   % an illumination distribution
estimate=dip_image(rand(sX,sY));
if (1) %calculate with the mask
    global myillu_mask;
    myillu_mask=cell(1,NumIm-1);
    tmp=newim(estimate);
    mykernel=rr(size(estimate))<2;
    mypos=floor(size(estimate)/2);
    tmp(mypos(1),mypos(2))=1;
    tmp(1,1)=1;
    tmp(10,10)=1;
    mymask=convolve(tmp,mykernel) > 0.5; %just one point in the middle
    %     myillu_mask=repmat(mymask,[1 1 NumIm]);
    for i=1:NumIm-1
        %     myillu_mask{i}=mymask;
        myillu_mask{i}=fft2rft(mymask);
    end
end
estimateRe=dip_image(rand([size(myillu_mask{1},2) size(myillu_mask{1},1)])); %real part
estimateIm=dip_image(rand([size(myillu_mask{1},2) size(myillu_mask{1},1)])); %imaginary part

%%   To not rerun the random generators

%disableCuda();
h=obj*0;h(2,2)=1;h=gaussf(h);
img=sqrt(prod(size(obj)))*real(ift(ft(obj) .* ft(h)));
%oimg=convolve(obj,h);

global myim;myim={};
for n=1:NumIm
    myim{n}=noise(img,'poisson'); 
end
global otfrep;otfrep={};otfrep{1}=rft(h);
global DeconvMethod;DeconvMethod='LeastSqr';
%global DeconvMethod;DeconvMethod='Poisson';
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
RegularisationParameters=ParseRegularisation(MyIlluReg,ToReg);
global ToEstimate;ToEstimate=1;   % 0 is object (with or without known illu), 1 is illu, 2 is OTF
AssignFunctions(RegularisationParameters,2)  % here: 0 is object, 1 is object with known illu, 2 is illum, 3 is OTF
%myVec=double(reshape(estimate,prod(size(estimate)))); % object estimate
% global OTFmask;
% OTFmask=PupilInterpolators.Mask;

% mysize=size(PupilInterpolators.indexList2D,2);
% mysize=[10 10];
% VecOTF=(rand(2*mysize,1)-0.5)*2;  % for real and imaginary part
% myVec=transpose(VecOTF);
% myVec=double(reshape(estimate(myillu_mask{1}), prod(size(estimate(myillu_mask{1}))))); % object estimate
myVec=[];
for n=1:numel(myillu_mask)
    myVecRe=double(reshape(estimateRe(myillu_mask{1}), prod(size(estimateRe(myillu_mask{1}))))); % object estimate
    myVecIm=double(reshape(estimateIm(myillu_mask{1}), prod(size(estimateIm(myillu_mask{1}))))); % object estimate
% myVec=complex(myVecRe,myVecIm);
    myNewVec=cat(2,myVecRe,myVecIm);
    myVec=cat(2,myVec,myNewVec);
end
% Old version:
% myVec=double(reshape(estimate,prod(size(estimate))));


global myillu_sumcond;
myillu_sumcond={NumIm};
% myillu_sumcond={};

global my_sumcond
my_sumcond={NumIm};

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
DeconvMask=cuda(DeconvMask);
aRecon=cuda(aRecon);

[err,grad]=GenericErrorAndDeriv(myVec);
end
%%  Now test the gradient in each dimension by numerically calculating the secans
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
% global savedATF;
% global myillu;
global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF
global mytmp;
mytmp={};
% AssignToGlobal(ConvertInputToModel(grad)); % changes otfrep{1} and SavedATF
% grad=savedATF;
AssignToGlobal(ConvertInputToModel(grad));
gradimg=cat(3,mytmp{:});
cat(3,mytmp{1},rft(myillu{1}))

AssignToGlobal(ConvertInputToModel(mygrad));
mygradimg=cat(3,mytmp{:});


if(0)
    grad= reshape(dip_image(grad','single'),size(img)); %bug if rft and mask
    % else
    % grad2=extract(dip_image(grad','single'),[prod(size(img))]); %because of the mask, length is too short
    % grad2=reshape(dip_image(grad2','single'),[ceil(sqrt(prod(size(grad2)))) ceil(sqrt(prod(size(grad2))))]);
    % % grad2= reshape(dip_image(grad','single'),[ceil(sqrt(prod(size(grad)))) ceil(sqrt(prod(size(grad))))]); %didnt work
    % % grad3=rft2fft(grad2); %should now be the same size as mygrad
    % grad=grad2;%to avoid renaming everywhere
    grad2=extract(dip_image(grad','single'),100); %because of the mask, length is too short
    grad3=reshape(dip_image(grad2','single'),[10 10]);
    grad=grad3;
    
    mygrad2=extract(dip_image(mygrad','single'),100); %because of the mask, length is too short
    mygrad3=reshape(dip_image(mygrad2','single'),[10 10]);
    mygrad=mygrad3;
end

% AssignToGlobal(ConvertInputToModel(mygrad)); % changes otfrep{1} and SavedATF
% mygrad=savedATF;
% mygrad=reshape(dip_image(mygrad','single'),size(img));


if isa(gradimg,'cuda')
    mygrad=cuda(mygradimg);
end

todisplay=cat(4,gradimg,mygradimg)

relerror = (mygrad - transpose(grad)) ./ mean(mean(abs(grad)));
fprintf('Max relative error :%g, \n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges

relerror = (mygradimg - gradimg) ./ mean(abs(gradimg),gradimg~=0);
fprintf('Max relative error :%g, \n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges
fprintf('Max relative center error :%g\n',max(abs(relerror(1:end-1,1:end-1,:))))   % Problems are caused by the hessian operator on finite arrays at the edges
