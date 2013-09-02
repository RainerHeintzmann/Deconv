% MyIdivErrorAndDeriv(aRecon) : Error measure for deconvolution algorithm
% this function interpretes the 4th dimension as a multi-view deconvolution
% each element of the 4th dimension having a corresponding PSF

%***************************************************************************
%   Copyright (C) 2008-2009 by Rainer Heintzmann                          *
%   heintzmann@gmail.com                                                  *
%                                                                         *
%   This program is free software; you can redistribute it and/or modify  *
%   it under the terms of the GNU General Public License as published by  *
%   the Free Software Foundation; Version 2 of the License.               *
%                                                                         *
%   This program is distributed in the hope that it will be useful,       *
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
%   GNU General Public License for more details.                          *
%                                                                         *
%   You should have received a copy of the GNU General Public License     *
%   along with this program; if not, write to the                         *
%   Free Software Foundation, Inc.,                                       *
%   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
%**************************************************************************
%

function [err,thegrad]=GenericErrorAndDeriv(myinput)
global myim;
global myillu;  % only of there is an illumination pattern present, will it be used
global otfrep;
global DeconvMethod;
global RegularisationParameters;  % This is a matrix with all possible regularisation lambdas (and other parameters)
global BetaVals;  % These are the scaling factors between pixel coordinates. This is proportional to the pixel width, but should be normalized 
global DeconvMask;  % only data in this mask will be evaluated
global DeconvVariances;  % only data in this mask will be evaluated
global ToEstimate;   % This flag controls what to estimate in this iteration step. 0: sample density, 1: illumination intensity, 2: both, 3: psf
global aRecon;   % Here the sample is stored, if the estimation is illumination or psf.
global NormFac;  % normalisation factor
global myillu_sumcond;   % denotes the positions in myillu for which each the sum condition (down to previous mention) are fullfilled.
global ForcePos;
global ComplexPSF; % If this flag is set, the PSF is complex and Fourier-transforms are done fully complex
global IntensityData; % If this flag is active, the data is interpreted as the abs square of the image amplitude (after convolution with the coherent asf)
global ComplexObj; % If active, the object is assumed to be complex valued. 

%delta= 100; % Weight for the negativtiy penalty
%delta= 1000; % Weight for the negativtiy penalty
%delta= 1; % Weight for the negativtiy penalty

DMask=[];
if ~isempty(DeconvMask)
    DMask=DeconvMask;
end

% Recast the matlab data into the dip_image datastructure
%aRecon=dip_image(aRecon);
norm3D = sqrt(prod(size(myim{1})));  % aRecon is not defined yet.
DataSize=size(myim{1});
DataLength=prod(DataSize);
%if length(DataSize) < 3
%    DataSize(3) = 1;
%end
if ForcePos && (isempty(ToEstimate) || ToEstimate==0)
    savedInput=myinput;
    myinput=abssqr(myinput); % converts the auxilary function back to the all positive object
end

if isempty(ToEstimate) || ToEstimate==0
    % aRecon=reshape(dip_image(myinput','single'),DataSize);  % The reconstruction can be up to 3D, whereas the data might be 4D
    aRecon=convertVecToObj(myinput,DataSize);  % converts the matlab vector into a dipimage object
    thegrad=newim(aRecon);  % clears the gradient. Defines it first as a dipimage. Later it is converted back to a double vector
elseif ToEstimate==1
    convertVecToIllu(myinput);  % will set the myillu{n} including the last one which is estimated from the previous ones
    asize=size(myillu{1});
    if length(asize) < 3
        asize(length(asize)+1:3)=1;
    end
    thegrad=newim([asize (length(myillu)-length(myillu_sumcond))]);  % clears the gradient. Defines it first as a dipimage. Later it is converted back to a double vector
else
    error('Other estimation methos not implemented yet.');
end
err=0;           % clears the errorsum
clear myinput;  

dip_setboundary('periodic')  % Establishes Periodic Boudary conditions.

currentSumCondIdx=1;   % goes through each of the sum conditions
prevSumCondGradIdx=0;

for viewNum = 1:length(myim)    % This iterates over all the different illumination patterns
myReg=0;
myRegGrad=0;
    
myImg=myim{viewNum};  % This does not cost any time or memory
myOtf=otfrep{1+mod(viewNum-1,length(otfrep))};  % This does not cost any time or memory. If only one otf is there it will always be used
if ~isempty(myillu)
    myIllum=myillu{1+mod(viewNum-1,length(myillu))};
else
    myIllum=1;
end

if ~isempty(myillu)
    if length(myillu) == 1  % this means the fft needs to be done only once
        if viewNum==1  % needs to be done only once in this case
            if ComplexPSF
                ftRecons=ft(aRecon .* myIllum);
            else
                ftRecons=rft(aRecon .* myIllum);
            end
        end
    else  % do a forward convolution for every instance of illumination distribution myillu
        if ComplexPSF
           ftRecons=ft(aRecon .* myIllum);
        else
           ftRecons=rft(aRecon .* myIllum);
        end
    end
    %if length(myillu) > 1
    %    error('For the moment only a single illumination distribution is allowed for each call to GenericErrorAndDeriv');
    %end
else  % no illumination is used
    % Forward model
    if viewNum==1  % needs to be done only once in this case
        if ComplexPSF
            ftRecons=ft(aRecon);
        else
            ftRecons=rft(aRecon);
        end
    end
end


%if isa(aRecon,'cuda')
if ComplexPSF
    Recons=norm3D*ift(ftRecons .* myOtf);  % convolve with corresponding PSF, still in Fourier space
else
    Recons=norm3D*rift(ftRecons .* myOtf);  % convolve with corresponding PSF, still in Fourier space
end
    
if IntensityData
    ReconsSaved=Recons;
    Recons=abssqr(Recons);  % After convolution the intensity is now calculated to be compared to the intensity data
end

% Residual: (here Cezar's i-divergence)
switch DeconvMethod
    case 'Poisson'
        % The first one below is numerically very unstable
        % myError =sum(Recons-myImg .* log(Recons));  % fast version of Czesar's i-divergence omitting constants
        eps=1e-7;  % To avoid to take log(0) in the Sterling approximation. If eps is too small, the ratio can generate quite high values in the ratio
        Recons(Recons<eps)=eps;
        ratio = myImg ./ Recons;
        epsR=1e-7;
        ratio(ratio<epsR)=epsR;
        %ratio(:,:,23)
        if ~isempty(DMask)
            ratio(~DMask)=1;  % This assumes always perfect agreement in this area
        end
        
        %myError =sum((Recons-myImg)+myImg .* log(ratio));  % fast version of Czesar's i-divergence omitting constants
        myError =sum((Recons-myImg)+myImg .* log(ratio),DMask);  % fast version of Czesar's i-divergence omitting constants
        %myError = 0;
        %myError =sum(Recons+myImg .* log(ratio));  % fast version of Czesar's i-divergence omitting constants
        %myError =sum(Recons-myImg .* log(Recons));  % fast version of Czesar's i-divergence omitting constants
        
        residuum=1-ratio;
        clear ratio
    case 'LeastSqr'
        residuum =( Recons-myImg);
        if ~isempty(DMask)
            residuum(~DMask)=0;
        end
        
        myError = real(sum(residuum .* conj(residuum),DMask));  % simple least squares, but accounting for complex numbers in residuum
        residuum=2*residuum; % to account for the square
    case 'WeightedLeastSqr'
        residuum =( Recons-myImg);
        if ~isempty(DMask)
            residuum(~DMask)=0;
        end
        
        if ~isempty(DeconvVariances)
            myError = real(sum(residuum .* conj(residuum) ./ DeconvVariances,DMask));  % simple least squares, but accounting for complex numbers in residuum
            residuum=2*residuum ./ DeconvVariances;
        else
            error('If using the WeightedLeastSqr method, you need ot provide a non-empty variance array');
        end
    otherwise
        error('Unknown decvonvolution method');
end

if IntensityData
    residuum = 2*ReconsSaved.*residuum;  % This treats goes back from the intensity world to the amplitude world
    clear ReconsSaved;
end

if ComplexPSF  % For all types of norms it is always the convolution with the conjugate OTF at the end
    myGrad = norm3D*ift(ft(residuum) .* conj(myOtf));
else
    myGrad = norm3D*rift(rft(residuum) .* conj(myOtf));
end

if ~ComplexObj
    myGrad=real(myGrad); % Is this correct? Simply discard the imaginary part here for objects which are known to be real?
end

clear myImg;
clear residuum;
clear myOtf;
clear Recons;

% The second part of the gradient calculation (multiplication with aRecon or illu{viewNum} is done at the end of the loop

% The terms below are only dependend on the reconstruction but not on the data. They need to be calculated onyl once for all views.

% Total error term (including regularisation):
% old version
% new version

if viewNum == 1 || ~isempty(ToEstimate) && ToEstimate==1
    if isempty(ToEstimate) || ToEstimate==0  
        toRegularize=aRecon;  % The regularisation needs to be applied to the object
    else
        toRegularize=myIllum;  % The regularisation needs to be applied to the illumination pattern
    end
end

myReg=0;
myRegGrad=0;
 
myLambda=RegularisationParameters(1,1); %case 'GS'  % Gradient roughness penalty: Gradient ^ 2
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGS(toRegularize,BetaVals);
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(2,1); %case 'AR'  % Arigovindan's roughness penalty
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeAR(toRegularize,BetaVals);
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end

myLambda=RegularisationParameters(3,1); %case 'TV'  Total variation roughness penalty (with eps)
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeTV(toRegularize,BetaVals,RegularisationParameters(3,2));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end

myLambda=RegularisationParameters(4,1); %case 'NegSqr'  Penalty for negative values
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeNegSqr(toRegularize);
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
    % fprintf('Neg Penalty: %g\n',delta*sum(aRecon.^2.*(aRecon<0))*lambdaPenalty /prod(size(aRecon)));
end
myLambda=RegularisationParameters(5,1); %case 'GR'  % Good's roughness penalty: Gradient ^ 2/f
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGR(toRegularize,BetaVals);
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end

clear ToRegularise;

if ~isempty(ToEstimate) && ToEstimate==1    % estimate the illumination
    myGrad= myGrad .*  aRecon;    % multiplication with the sample density (rho) for estimating the gradient of myillu
    agradIdx=viewNum-1-(currentSumCondIdx-1);
    if viewNum ~= myillu_sumcond{currentSumCondIdx}
        thegrad(:,:,:,agradIdx)=myGrad + myRegGrad;
    else  % The last residuum has to be subtracted from each of the other residuals, see eq. S14 and S4 in supplementary methods of DOI: 10.1038/NPHOTON.2012.83
        thegrad(:,:,:,prevSumCondGradIdx:agradIdx-1)=thegrad(:,:,:,prevSumCondGradIdx:agradIdx-1)-repmat(myGrad+myRegGrad,[1 1 1 agradIdx-prevSumCondGradIdx]);
        prevSumCondGradIdx=agradIdx;
        currentSumCondIdx = currentSumCondIdx +1;
    end
    
    err=err+myError + myReg;
else
    if ~isempty(myillu)   % account for the illumination pattern in the object iteration (if present)
        myGrad=myGrad .*  myIllum;
    end
    if viewNum == 1       % object penalty needs to be accounted for only once
        thegrad=thegrad+ myGrad +  myRegGrad; 
        err=err+myError + myReg;    % was cleared before the for-loop
    else
        thegrad=thegrad+myGrad;   % was cleared before the for-loop
        err=err+myError;    % was cleared before the for-loop
    end
end

end  % viewNum
clear ftRecons;
clear DMask;
clear myGrad;

err = double(NormFac*err);

%dipshow(7,grad(:,:,7))
% Line below also takes care of possible crunch operation (ToEstimate is 1 and myillu_mask exists
thegrad=convertGradToVec(NormFac*thegrad);    % converts the dip_image back to a linear matlab vector

if ForcePos && (isempty(ToEstimate) || ToEstimate==0)
    thegrad = 2*savedInput.*thegrad;  % To account for the fact that the auxilary function is what is iterated and the object estimate is the square of it
end

%fprintf('Val: %g, Penalty %g, Gradien Norm: %g Penalty %g\n',err, lambdaPenalty*myReg,norm(grad),norm(lambdaPenalty * double(myRegGrad(:))));
%fprintf('Val: %g, Penalty %g\n',err, lambdaPenalty*myReg);

%fprintf('X: %g, Total Grad: %g\n', sum(abs(aRecon)),sum(abs(grad)))

%grad=single_force(grad);
