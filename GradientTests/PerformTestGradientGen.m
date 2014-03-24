%% Generic test gradient function
% In should be a list of all the types that should be tested with options (for example obj, illumination, with or without mask, with or without ForcePos...)
% What about TestObjGradientIllu???
% Out should be the error for this case
%
% PARAMETERS are (possibilities):
% myType: type of variable to test (obj, illu, reg, otf). Default obj
% MyReg: type of regularization (same syntax and options as for GenericDeconvolution, exept for empty regularization, in which case write 0). Default 0
% useMask: only for type=illu. (1, 0). If 1, uses a pre-defined mask with 3 areas. If 0, myillu_mask=[]; Default 0
% useCuda (1, 0). Default 0
% dbg_display: Displays the comparison between nominal and calculated gradient (1, 0). Default 1
% NumIm: number of images. Default 3
% sX: size along X. Default 102
% sY: size along Y. Default 100

function out=PerformTestGradientGen(varargin)

%% Input arguments

% Define defaults
options = struct('myType','obj','MyReg',0,'useMask',0,'useCuda',0,...
    'dbg_display',1,'NumIm',3,'sX',102,'sY', 100);

% read the acceptable names
optionNames = fieldnames(options);

% count arguments
nArgs = length(varargin);
if round(nArgs/2)~=nArgs/2
   error('PerformTestGradientGen needs paramName/paramValue pairs')
end

for pair = reshape(varargin,2,[]) % pair is {paramName;paramValue}
   inpName = pair{1}; 

   if any(strmatch(inpName,optionNames)) % looks for matching parameter names
      options.(inpName) = pair{2}; % overwrite parameters
   else
      error('%s is not a recognized parameter name',inpName)
   end
end

for i=1:length(optionNames)
    str=[optionNames{i} '=options.' optionNames{i} ';']; % for example obj=options.obj
    eval(str); % executes obj=options.obj;
end

%% Common code

disableCuda();

% Create images that are common to all cases
% Not sure about the random seeds
rng(1); % initialize the raqndom generator with the same seed always
obj=dip_image(rand(sX,sY));
rng(2); % initialize the raqndom generator with the same seed always
objestimate=10*dip_image(rand(sX,sY));  % current (old) reconstruction estimate
rng(3); % initialize the raqndom generator with the same seed always
illu=dip_image(rand(sX,sY,NumIm));   % an illumination distribution
rng(4); % initialize the raqndom generator with the same seed always
estimate=dip_image(rand(sX,sY));

h=obj*0;h(2,2)=1;h=gaussf(h);
img=sqrt(prod(size(obj)))*real(ift(ft(obj) .* ft(h)));
%oimg=convolve(obj,h);

global myim;myim={};
for n=1:NumIm
%     rng(n+6); % initialize the raqndom generator with the same seed always
%     myim{n}=noise(img,'poisson'); 
    myim{n}=img; %to allow comparison 
end
global otfrep;otfrep={};otfrep{1}=rft(h);
global DeconvMethod;DeconvMethod='LeastSqr';
%global DeconvMethod;DeconvMethod='Poisson';
global DeconvVariance;
DeconvVariances=[];
global BetaVals;
BetaVals=[1 1 1];
global NegPenalty;NegPenalty='NONE'; %what is this?
global RegularisationMethod;RegularisationMethod='NONE'; %what is this?
%global DeconvMask;DeconvMask=xx(sX,sY)>0;  % only data in this mask will be evaluated
global DeconvMask;DeconvMask=[];  % only data in this mask will be evaluated
%global myillu;myillu{1}=illu;  % only data in this mask will be evaluated
global myillu;myillu={}; % myillu{1}=illu(:,:,0);myillu{2}=illu(:,:,1); 
global NormFac;NormFac=1.0;   % Normalisation factor
global aRecon;aRecon=objestimate;

% global PupilInterpolators;
% PupilInterpolators={};
global RegularisationParameters;

global my_sumcond
my_sumcond={NumIm};

if ~iscell(MyReg) && MyReg==0
    MyReg={}; %Because it not allowed as input in a varargin function.
end

%% Fill in options

switch myType
    
    
    case 'obj' %Object gradient estimation
        if useMask==1
            fprintf('useMask option set to 1 whereas this is not illumination estimation case! Ignoring this field and proceeding');
        end
        
        
    case 'illu' %Illu gradient estimation
        global myillu_mask;
        myillu_mask={};
        global myillu_sumcond;
        myillu_sumcond={NumIm};
        
        if useMask==1
            fprintf('An illu mask is used\n');
            myillu_mask=cell(1,NumIm-1);
            tmp=newim(estimate);
            mykernel=rr(size(estimate))<2;
            mypos=floor(size(estimate)/2);
            tmp(mypos(1),mypos(2))=1;
            tmp(1,1)=1;
            tmp(10,10)=1;
            mymask=convolve(tmp,mykernel) > 0.5; %just one point in the middle
            %     myillu_mask=repmat(mymask,[1 1 NumIm]);
            for n=1:NumIm-1
                %     myillu_mask{n}=mymask;
                myillu_mask{n}=fft2rft(mymask);
            end
            
            global mytmp;
            mytmp={};  % Used for debugging.
            
            rng(5); % initialize the raqndom generator with the same seed always
            estimateRe=dip_image(rand([size(myillu_mask{1},2) size(myillu_mask{1},1)])); %real part
            rng(6); % initialize the raqndom generator with the same seed always
            estimateIm=dip_image(rand([size(myillu_mask{1},2) size(myillu_mask{1},1)])); %imaginary part
            
            myVec=[];
            for n=1:numel(myillu_mask)
                myVecRe=double(reshape(estimateRe(myillu_mask{1}), prod(size(estimateRe(myillu_mask{1}))))); % object estimate
                myVecIm=double(reshape(estimateIm(myillu_mask{1}), prod(size(estimateIm(myillu_mask{1}))))); % object estimate
            % myVec=complex(myVecRe,myVecIm);
                myNewVec=cat(2,myVecRe,myVecIm);
                myVec=cat(2,myVec,myNewVec);
            end
        end %end case of illu mask
        
        ToReg=1;  % 0 is object, 1 means illu, 2 means otf
        RegularisationParameters=ParseRegularisation(MyReg,ToReg);
        global ToEstimate;ToEstimate=1;   % 0 is object (with or without known illu), 1 is illu, 2 is OTF
        AssignFunctions(RegularisationParameters,2)  % here: 0 is object, 1 is object with known illu, 2 is illum, 3 is OTF
    
    
    case 'reg' %Regularizer gradient estimation
        if useMask==1
            fprintf('useMask option set to 1 whereas this is not illumination estimation case! Ignoring this field and proceeding');
        end
        
        
    case 'otf' %For Rainer
        
end

%% Now calculate the gradient

disableCuda();  % Always run this part without cuda.
[err,grad]=GenericErrorAndDeriv(myVec);  % to get the analytical gradient solution
dipgrad=grad; diperr=err; dipVec=myVec;

%%   Do all the cuda conversions here to test the cuda variant of the algorithm
if (useCuda)
initCuda();
enableCuda()
% myVec=cuda(myVec);
myVec=ConditionalCudaConvert(myVec,useCuda);
% grad=cuda(grad);
grad=ConditionalCudaConvert(grad,useCuda);
% for n=1:numel(myillu)
%     myillu{n}=cuda(myillu{n});
%     myim{n}=cuda(myim{n});
% end
myillu=ConditionalCudaConvert(myillu,useCuda);
myim=ConditionalCudaConvert(myim,useCuda);

otfrep{1}=rft(fftshift(cuda(h)));
% DeconvMask=cuda(DeconvMask);
DeconvMask=ConditionalCudaConvert(DeconvMask,useCuda);
% for n=1:NumIm-1
%     myillu_mask{n}=cuda(myillu_mask{n});
% end
myillu_mask=ConditionalCudaConvert(myillu_mask,useCuda);
% aRecon=cuda(aRecon);
aRecon=ConditionalCudaConvert(aRecon,useCuda);

[err,grad]=GenericErrorAndDeriv(myVec);
relErr=(diperr-err)/err;
relGradErr=max(abs(dipgrad-grad))/mean(abs(grad));
VecErr=max(abs(dipVec-myVec));
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

% max(abs(transpose(dipgrad)-mygrad))/mean(abs(mygrad))
% max(abs(transpose(grad)-mygrad))/mean(abs(mygrad))


global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF

AssignToGlobal(ConvertInputToModel(grad));
gradimg=cat(3,mytmp{:});
% cat(3,mytmp{1},rft(myillu{1}))

AssignToGlobal(ConvertInputToModel(mygrad));
mygradimg=cat(3,mytmp{:});

useCuda=0; disableCuda();
if dbg_display==1
    global todisplay
    todisplay=cat(4,gradimg,mygradimg)
end

% relerror = (mygrad - transpose(grad)) ./ mean(mean(abs(grad)));
% fprintf('Max relative error :%g, \n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges

relerror = (mygradimg - gradimg) ./ mean(abs(gradimg),gradimg~=0);
% fprintf('Max relative error :%g, \n',max(abs(relerror)))   % Problems are caused by the hessian operator on finite arrays at the edges
% fprintf('Max relative center error :%g\n',max(abs(relerror(1:end-1,1:end-1,:))))   % Problems are caused by the hessian operator on finite arrays at the edges

out=max(abs(relerror));

end