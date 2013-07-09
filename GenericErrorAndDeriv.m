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
% This function is based on a talk by 
% Muthuvel  Arigovindan*,  Daniel Elnatan,  Jennifer  Fung,  Eric Branlund,  John W. Sedat,   and  David A.  Agard
% University of California at San Francisco
% suggesting a regolarisation term as given below

function [err,grad]=GenericErrorAndDeriv(myinput)
global myim;
global myillu;  % only of there is an illumination pattern present, will it be used
global otfrep;
global lambdaPenalty;
global DeconvMethod;
global RegularisationMethod;
global NegPenalty;
global BetaVals;  % These are the scaling factors between pixel coordinates. This is proportional to the pixel width, but should be normalized 
global DeconvMask;  % only data in this mask will be evaluated
global DeconvVariances;  % only data in this mask will be evaluated
global ToEstimate;   % This flag controls what to estimate in this iteration step. 0: sample density, 1: illumination intensity, 2: both, 3: psf
global aRecon;   % Here the sample is stored, if the estimation is illumination or psf.
global NormFac;  % normalisation factor
global myillu_sumcond;   % denotes the positions in myillu for which each the sum condition (down to previous mention) are fullfilled.

alpha= 1.0;  % Weighting between laplacian and noise penalty
%alpha= 10000.0;  % Weighting between laplacian and noise penalty
gamma = 1.0; % Weight for the intensity penalty under the logarithm
%delta= 100; % Weight for the negativtiy penalty
delta= 1000; % Weight for the negativtiy penalty
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

if isempty(ToEstimate) || ToEstimate==0
    % aRecon=reshape(dip_image(myinput','single'),DataSize);  % The reconstruction can be up to 3D, whereas the data might be 4D
    aRecon=convertVecToObj(myinput,DataSize);  % converts the matlab vector into a dipimage object
    grad=newim(aRecon);  % clears the gradient. Defines it first as a dipimage. Later it is converted back to a double vector
elseif ToEstimate==1
    convertVecToIllu(myinput);  % will set the myillu{n} including the last one which is estimated from the previous ones
    asize=size(myillu{1});
    if length(asize) < 3
        asize(length(asize)+1:3)=1;
    end
    grad=newim([asize (length(myillu)-length(myillu_sumcond))]);  % clears the gradient. Defines it first as a dipimage. Later it is converted back to a double vector
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
            %if isa(aRecon,'cuda')
                ftRecons=rft(aRecon .* myIllum);
            %else
            %    ftRecons=ft(aRecon .* myIllum);
            %end
        end
    else  % do a forward convolution for every instance of illumination distribution myillu
        %if isa(aRecon,'cuda')
            ftRecons=rft(aRecon .* myIllum);
        %else
        %    ftRecons=ft(aRecon .* myIllum);
        %end
    end
    %if length(myillu) > 1
    %    error('For the moment only a single illumination distribution is allowed for each call to GenericErrorAndDeriv');
    %end
else  % no illumination is used
    % Forward model
    if viewNum==1  % needs to be done only once in this case
        %if isa(aRecon,'cuda')
            ftRecons=rft(aRecon);
        %else
        %    ftRecons=ft(aRecon);
        %end
    end
end


%if isa(aRecon,'cuda')
    Recons=norm3D*rift(ftRecons .* myOtf);  % convolve with corresponding PSF, still in Fourier space
%else
%    Recons=norm3D*real(ift(ftRecons .* myOtf));  % convolve with corresponding PSF, still in Fourier space
%end
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
        
        %if isa(aRecon,'cuda')
            myGrad = norm3D*rift(rft(1-ratio) .* conj(myOtf));
        %else
        %    myGrad = norm3D*real(ift(ft(1-ratio) .* conj(myOtf)));
        %end
        clear ratio
    case 'LeastSqr'
        residuum =( Recons-myImg);
        if ~isempty(DMask)
            residuum(~DMask)=0;
        end
        
        myError = real(sum(residuum .* conj(residuum),DMask));  % simple least squares, but accounting for complex numbers in residuum
        %if isa(aRecon,'cuda')
            myGrad = norm3D*2*rift(rft(conj(residuum)) .* conj(myOtf));
        %else
        %    myGrad = norm3D*2*real(ift(ft(conj(residuum)) .* conj(myOtf)));
        %end
    case 'WeightedLeastSqr'
        residuum =( Recons-myImg);
        if ~isempty(DMask)
            residuum(~DMask)=0;
        end
        
        if ~isempty(DeconvVariances)
            myError = real(sum(residuum .* conj(residuum) ./ DeconvVariances,DMask));  % simple least squares, but accounting for complex numbers in residuum
            
            %if isa(aRecon,'cuda')
                myGrad = norm3D*2*rift(rft(conj(residuum ./ DeconvVariances)) .* conj(myOtf));
            %else
            %    myGrad = norm3D*2*real(ift(ft(conj(residuum ./ DeconvVariances)) .* conj(myOtf)));
            %end
        else
            error('If using the WeightedLeastSqr method, you need ot provide a non-empty variance array');
        end
    otherwise
        error('Unknown decvonvolution method');
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

switch RegularisationMethod
    case 'NONE'
        myReg=0;
        myRegGrad=0;
    case 'GR'  % Good's roughness penalty: Abs(Gradient ^ 2)
        
        % aGrad=gradient(aRecon)
        % myHessian = hessian(aRecon);
        % myRegGrad = -2*real(myHessian{1,1}+myHessian{2,2}+2*myHessian{1,2}); 
        if (ndims(toRegularize) == 2) || (size(toRegularize,3) == 1)
            if (ndims(toRegularize) == 2)
                aGrad{1}=(circshift(toRegularize,[1 0])-circshift(toRegularize,[-1 0]))/(2*BetaVals(1));  % cyclic rotation
                aGrad{2}=(circshift(toRegularize,[0 1])-circshift(toRegularize,[0 -1]))/(2*BetaVals(2));  % cyclic rotation
                myRegGrad = (1/(2*BetaVals(1)*BetaVals(1)))*(2*toRegularize-circshift(toRegularize,[-2 0]) - circshift(toRegularize,[2 0])) + ...
                    (1/(2*BetaVals(2)*BetaVals(2)))*(2*toRegularize-circshift(toRegularize,[0 -2]) - circshift(toRegularize,[0 2]));
            else (size(toRegularize,3) == 1)
                aGrad{1}=(circshift(toRegularize,[1 0 0])-circshift(toRegularize,[-1 0 0]))/2;  % cyclic rotation
                aGrad{2}=(circshift(toRegularize,[0 1 0])-circshift(toRegularize,[0 -1 0]))/2;  % cyclic rotation
                myRegGrad = (1/(2*BetaVals(1)*BetaVals(1)))*(2*toRegularize-circshift(toRegularize,[-2 0 0]) - circshift(toRegularize,[2 0 0])) + ...
                    (1/(2*BetaVals(2)*BetaVals(2)))*(2*toRegularize-circshift(toRegularize,[0 -2 0]) - circshift(toRegularize,[0 2 0]));
            end
        
            % myReg = sum(BetaVals(1)*aGrad{1} .* aGrad{1} + BetaVals(2)*aGrad{2} .* aGrad{2});
            myReg = sum(aGrad{1} .* aGrad{1} + aGrad{2} .* aGrad{2});
        elseif ndims(toRegularize) == 3
            aGrad{1}=(circshift(toRegularize,[1 0 0])-circshift(toRegularize,[-1 0 0]))/(2*BetaVals(1));  % cyclic rotation
            aGrad{2}=(circshift(toRegularize,[0 1 0])-circshift(toRegularize,[0 -1 0]))/(2*BetaVals(2));  % cyclic rotation
            aGrad{3}=(circshift(toRegularize,[0 0 1])-circshift(toRegularize,[0 0 -1]))/(2*BetaVals(3));  % cyclic rotation
            myReg = sum(aGrad{1} .* aGrad{1} + aGrad{2} .* aGrad{2} + aGrad{3} .* aGrad{3});
            myRegGrad = (1/(2*BetaVals(1)*BetaVals(1)))*(2*toRegularize-circshift(toRegularize,[-2 0 0])-circshift(toRegularize,[2 0 0])) + ...
                        (1/(2*BetaVals(2)*BetaVals(2)))*(2*toRegularize-circshift(toRegularize,[0 -2 0])-circshift(toRegularize,[0 2 0])) + ...
                        (1/(2*BetaVals(3)*BetaVals(3)))*(2*toRegularize-circshift(toRegularize,[0 0 -2])-circshift(toRegularize,[0 0 2]));
        else % 1-D
            aGrad{1}=(circshift(toRegularize,1)-circshift(toRegularize,-1))/(2*BetaVals(1));  % cyclic rotation
            myReg = sum(aGrad{1} .* aGrad{1});
            myRegGrad = (1/(2*BetaVals(1)*BetaVals(1)))*(2*toRegularize-circshift(toRegularize,-2) - circshift(toRegularize,2)); 
        end            
        
    case 'AR'
        % H = hessian(aRecon);
        aReconSqr=toRegularize.*toRegularize;

        if ndims(toRegularize) == 2 || (size(toRegularize,3) == 1)
            s1=1/1.6;s2=0.1727;
            H11=dip_convolve1d(toRegularize,s1*[1 -2 1],0,1)/(BetaVals(1)*BetaVals(1)); % second X derivative
            H22=dip_convolve1d(toRegularize,s1*[1 -2 1],1,1)/(BetaVals(2)*BetaVals(2)); % second Y derivative
            H2a=dip_convolve1d(toRegularize,[-1 0 1],0,1); %
            H12=dip_convolve1d(H2a,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
            myProjHessianSqr= H11 .*H11 + H22.*H22 + 2*H12.*H12; 
            T1 = 1+ alpha*(gamma*aReconSqr + myProjHessianSqr);  % will be summed over as log values
            tmp11 = dip_convolve1d(H11 ./ T1,s1*[1 -2 1],0,1)/(BetaVals(1)*BetaVals(1));
            tmp12=dip_convolve1d(H12 ./ T1,[-1 0 1],0,1); %
            tmp12=2*dip_convolve1d(tmp12,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
            tmp22 = dip_convolve1d((H22 ./ T1),s1*[1 -2 1],1,1)/(BetaVals(2)*BetaVals(2));
            myHessian2 = tmp11 + tmp12 + tmp22;
            myRegGrad = 2*alpha*myHessian2 + gamma*2*alpha*toRegularize ./ T1;
        else
            s1=1/1.6;s2=0.1727;
            H11=dip_convolve1d(toRegularize,s1*[1 -2 1],0,1)/(BetaVals(1)*BetaVals(1)); % second X derivative
            H22=dip_convolve1d(toRegularize,s1*[1 -2 1],1,1)/(BetaVals(2)*BetaVals(2)); % second Y derivative
            H33=dip_convolve1d(toRegularize,s1*[1 -2 1],2,1)/(BetaVals(3)*BetaVals(3)); % second Y derivative
            H2a=dip_convolve1d(toRegularize,[-1 0 1],0,1); %
            H12=dip_convolve1d(H2a,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
            H2b=dip_convolve1d(toRegularize,[-1 0 1],0,1); %
            H13=dip_convolve1d(H2b,s2*[-1 0 1],2,1)/(BetaVals(1)*BetaVals(3)); % mixed derivative XZ
            H2c=dip_convolve1d(toRegularize,[-1 0 1],1,1); %
            H23=dip_convolve1d(H2c,s2*[-1 0 1],2,1)/(BetaVals(2)*BetaVals(3)); % mixed derivative YZ
            myProjHessianSqr= H11 .*H11 + H22.*H22 +H33.*H33 + 2*H12.*H12+ 2*H13.*H13+ 2*H23.*H23;  % Does it need the weights of the filters?
            T1 = 1+ alpha*(gamma*aReconSqr + myProjHessianSqr);
            tmp11 = dip_convolve1d(H11 ./ T1,s1*[1 -2 1],0,1)/(BetaVals(1)*BetaVals(1));
            tmp22 = dip_convolve1d(H22 ./ T1,s1*[1 -2 1],1,1)/(BetaVals(2)*BetaVals(2));
            tmp33 = dip_convolve1d(H33 ./ T1,s1*[1 -2 1],1,1)/(BetaVals(3)*BetaVals(3));
            tmp12=dip_convolve1d(H12 ./ T1,[-1 0 1],0,1); %
            tmp12=2*dip_convolve1d(tmp12,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
            tmp13=dip_convolve1d(H13 ./ T1,[-1 0 1],0,1); %
            tmp13=2*dip_convolve1d(tmp13,s2*[-1 0 1],2,1)/(BetaVals(1)*BetaVals(3)); % mixed derivative XZ
            tmp23=dip_convolve1d(H23 ./ T1,[-1 0 1],1,1); %
            tmp23=2*dip_convolve1d(tmp23,s2*[-1 0 1],2,1)/(BetaVals(2)*BetaVals(3)); % mixed derivative YZ

            myHessian2 = tmp11 + tmp12 + tmp13 + tmp22+ tmp23+tmp33 ;
            myRegGrad = 2*alpha*myHessian2 + gamma*2*alpha*toRegularize ./ T1;
        end
        clear aReconSqr;
        % An alternative would have been to use the theorem:
        % (d/df)(d/dx)f(x) = f''(x) / f'(x)
        
        myReg = sum(log(T1));
    case 'TV'
        % The total variation code below is based on a simple 2-point finite difference
        % calculation using cicular shifting (like rsl, rsr)
        % The derivative with respect to the function value is calculated
        % as the exact value for the discrete version of the finite
        % differences. This way the gradient should be correct.
        % It was checked using the TestGradients function.
        % aGrad=gradient(aRecon);
        % betaX=1.0;betaY=1.0;betaZ=1.0; % In the future this can be anisotropic
        %
        % The Regularisation modification with epsR was introduced, according to 
        % Ferreol Soulez et al. "Blind deconvolution of 3D data in wide field fluorescence microscopy
        % 
        epsC=1e-10;
        epsR=1e2*1e2;  % Modifies the TV norm to behave better (e.g. considering digitalisation effects)
        if (ndims(toRegularize) == 2) || (size(toRegularize,3) == 1)
            if (ndims(toRegularize) == 2)
                aGrad{1}=(circshift(toRegularize,[1 0])-circshift(toRegularize,[-1 0]))/(2*BetaVals(1));  % cyclic rotation
                aGrad{2}=(circshift(toRegularize,[0 1])-circshift(toRegularize,[0 -1]))/(2*BetaVals(2));  % cyclic rotation
                mySqrt = sqrt(aGrad{1}*aGrad{1}+aGrad{2}*aGrad{2}+epsR);
                mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
                myRegGrad = (1/(4*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,[-2 0]))./circshift(mySqrt,[-1 0]) - (circshift(toRegularize,[2 0])-toRegularize)./circshift(mySqrt,[1 0])) + ...
                        (1/(4*BetaVals(2)*BetaVals(2)))*((toRegularize-circshift(toRegularize,[0 -2]))./circshift(mySqrt,[0 -1]) - (circshift(toRegularize,[0 2])-toRegularize)./circshift(mySqrt,[0 1]));
            else (size(toRegularize,3) == 1)
                aGrad{1}=(circshift(toRegularize,[1 0 0])-circshift(toRegularize,[-1 0 0]))/(2*BetaVals(1));  % cyclic rotation
                aGrad{2}=(circshift(toRegularize,[0 1 0])-circshift(toRegularize,[0 -1 0]))/(2*BetaVals(2));  % cyclic rotation
                mySqrt = sqrt(aGrad{1}*aGrad{1}+aGrad{2}*aGrad{2}+epsR);
                mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
                myRegGrad = (1/(4*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,[-2 0 0]))./circshift(mySqrt,[-1 0 0]) - (circshift(toRegularize,[2 0 0])-toRegularize)./circshift(mySqrt,[1 0 0])) + ...
                        (1/(4*BetaVals(2)*BetaVals(2)))*((toRegularize-circshift(toRegularize,[0 -2 0]))./circshift(mySqrt,[0 -1 0]) - (circshift(toRegularize,[0 2 0])-toRegularize)./circshift(mySqrt,[0 1 0]));
            end
        elseif ndims(toRegularize) == 3
            aGrad{1}=(circshift(toRegularize,[1 0 0])-circshift(toRegularize,[-1 0 0]))/(2*BetaVals(1));  % cyclic rotation
            aGrad{2}=(circshift(toRegularize,[0 1 0])-circshift(toRegularize,[0 -1 0]))/(2*BetaVals(2));  % cyclic rotation
            aGrad{3}=(circshift(toRegularize,[0 0 1])-circshift(toRegularize,[0 0 -1]))/(2*BetaVals(3));  % cyclic rotation
            mySqrt = sqrt(aGrad{1}*aGrad{1}+aGrad{2}*aGrad{2}+aGrad{3}*aGrad{3}+epsR);
            mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
            myRegGrad = (1/(4*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,[-2 0 0]))./circshift(mySqrt,[-1 0 0]) - (circshift(toRegularize,[2 0 0])-toRegularize)./circshift(mySqrt,[1 0 0])) + ...
                        (1/(4*BetaVals(2)*BetaVals(2)))*((toRegularize-circshift(toRegularize,[0 -2 0]))./circshift(mySqrt,[0 -1 0]) - (circshift(toRegularize,[0 2 0])-toRegularize)./circshift(mySqrt,[0 1 0])) + ...
                        (1/(4*BetaVals(3)*BetaVals(3)))*((toRegularize-circshift(toRegularize,[0 0 -2]))./circshift(mySqrt,[0 0 -1]) - (circshift(toRegularize,[0 0 2])-toRegularize)./circshift(mySqrt,[0 0 1]));
        else % 1-D
            aGrad{1}=(circshift(toRegularize,1)-circshift(toRegularize,-1))/(2*BetaVals(1));  % cyclic rotation
            mySqrt = sqrt(aGrad{1}*aGrad{1}+epsR);
            mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
            myRegGrad = (1/(4*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,-2))./circshift(mySqrt,-1) - (circshift(toRegularize,2)-toRegularize)./circshift(mySqrt,1)); 
        end            
        myReg = sum(mySqrt);
    otherwise
        error('Unknown regularisation method: Options are: NONE, TV, AR');
end

switch NegPenalty
    case 'NONE'
    case 'NegSqr'
            myReg = myReg+delta*sum(toRegularize.^2.*(toRegularize<0));
            myRegGrad = myRegGrad+2*delta*toRegularize .* (toRegularize<0);
        % fprintf('Neg Penalty: %g\n',delta*sum(aRecon.^2.*(aRecon<0))*lambdaPenalty /prod(size(aRecon)));
    otherwise
        error('Unknown negative penalty method: Options are: NONE, NegSqr');
end

clear ToRegularise;

if ~isempty(ToEstimate) && ToEstimate==1    % estimate the illumination
    myGrad= myGrad .*  aRecon;    % multiplication with the sample density (rho) for estimating the gradient of myillu
    agradIdx=viewNum-1-(currentSumCondIdx-1);
    if viewNum ~= myillu_sumcond{currentSumCondIdx}
        grad(:,:,:,agradIdx)=myGrad + lambdaPenalty * myRegGrad;
    else  % The last residuum has to be subtracted from each of the other residuals, see eq. S14 and S4 in supplementary methods of DOI: 10.1038/NPHOTON.2012.83
        grad(:,:,:,prevSumCondGradIdx:agradIdx-1)=grad(:,:,:,prevSumCondGradIdx:agradIdx-1)-repmat(myGrad+lambdaPenalty * myRegGrad,[1 1 1 agradIdx-prevSumCondGradIdx]);
        prevSumCondGradIdx=agradIdx;
        currentSumCondIdx = currentSumCondIdx +1;
    end
    
    err=err+myError + lambdaPenalty * myReg;
else
    if ~isempty(myillu)   % account for the illumination pattern in the object iteration (if present)
        myGrad=myGrad .*  myIllum;
    end
    if viewNum == 1       % object penalty needs to be accounted for only once
        grad=grad+ myGrad + lambdaPenalty * myRegGrad; 
        err=err+myError + lambdaPenalty * myReg;    % was cleared before the for-loop
    else
        grad=grad+myGrad;   % was cleared before the for-loop
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
grad=convertGradToVec(NormFac*grad);    % converts the dip_image back to a linear matlab vector

%fprintf('Val: %g, Penalty %g, Gradien Norm: %g Penalty %g\n',err, lambdaPenalty*myReg,norm(grad),norm(lambdaPenalty * double(myRegGrad(:))));
%fprintf('Val: %g, Penalty %g\n',err, lambdaPenalty*myReg);

%fprintf('X: %g, Total Grad: %g\n', sum(abs(aRecon)),sum(abs(grad)))

%grad=single_force(grad);
