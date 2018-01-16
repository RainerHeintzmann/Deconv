% [myReg,myRegGrad]=RegularizeKevran(toRegularize,BetaVals,epsR) computes Total variation regularisation
% epsR [0..1] shifts the weight between Data and gradient regularization
% penalty = sqrt(epsR*|grad(f)|^2+(1-epsR)*f.^2)
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient
% The code below is based on a simple 2-point finite difference
% calculation using circular shifting (like rsl, rsr)
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

function [myReg,myRegGrad]=RegularizeKevran(toRegularize,BetaVals,epsR)
epsC=1e-10;
if nargin < 3
    epsR=0.5;  % Modifies the TV norm to behave better (e.g. considering digitalisation effects)
end
toRegularizeSqr=toRegularize.^2;
if (ndims(toRegularize) == 2) || (size(toRegularize,3) == 1)
    if (ndims(toRegularize) == 2)
        aGradL{1}=(toRegularize - circshift(toRegularize,[1 0]))/BetaVals(1);  % cyclic rotation
        aGradL{2}=(toRegularize - circshift(toRegularize,[0 1]))/BetaVals(2);  % cyclic rotation
        aGradR{1}=(circshift(toRegularize,[-1 0]) - toRegularize)/BetaVals(1);
        aGradR{2}=(circshift(toRegularize,[0 -1]) - toRegularize)/BetaVals(2);
        mySqrtL = sqrt(epsR*(abssqr(aGradL{1}) + abssqr(aGradL{2})) + (1.0-epsR)*toRegularizeSqr);
        mySqrtR = sqrt(epsR*(abssqr(aGradR{1}) + abssqr(aGradR{2})) + (1.0-epsR)*toRegularizeSqr);
        mySqrt = mySqrtL + mySqrtR;
        mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
        % gradient is tested
        myRegGrad = (epsR*(aGradL{1}/BetaVals(1) + aGradL{2}/BetaVals(2))+ (1.0-epsR)*2*toRegularize)./mySqrtL - (epsR*(aGradR{1}/BetaVals(1) + aGradR{2}/BetaVals(2))+ (1.0-epsR)*2*toRegularize)./mySqrtR - ...
                    aGradR{1}/BetaVals(1)./circshift(mySqrtL, [-1 0]) - aGradR{2}/BetaVals(2)./circshift(mySqrtL, [0 -1]) + ...
                    aGradL{1}/BetaVals(1)./circshift(mySqrtR, [1 0]) + aGradL{2}/BetaVals(2)./circshift(mySqrtR, [0 1]);
    else (size(toRegularize,3) == 1)
        aGradL{1}=(toRegularize - circshift(toRegularize,[1 0 0]))/BetaVals(1);  % cyclic rotation
        aGradL{2}=(toRegularize - circshift(toRegularize,[0 1 0]))/BetaVals(2);  % cyclic rotation
        aGradR{1}=(circshift(toRegularize,[-1 0 0]) - toRegularize)/BetaVals(1);
        aGradR{2}=(circshift(toRegularize,[0 -1 0]) - toRegularize)/BetaVals(2);
        mySqrtL = sqrt(epsR*(abssqr(aGradL{1})+abssqr(aGradL{2}))+ (1.0-epsR)*toRegularizeSqr);
        mySqrtR = sqrt(epsR*(abssqr(aGradR{1})+abssqr(aGradR{2}))+ (1.0-epsR)*toRegularizeSqr);
        mySqrt =  mySqrtL + mySqrtR;
        mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
        myRegGrad = (epsR*(aGradL{1}/BetaVals(1) + aGradL{2}/BetaVals(2))+ (1.0-epsR)*2*toRegularize)./mySqrtL - (epsR*(aGradR{1}/BetaVals(1) + aGradR{2}/BetaVals(2))+ (1.0-epsR)*2*toRegularize)./mySqrtR - ...
                    aGradR{1}/BetaVals(1)./circshift(mySqrtL, [-1 0 0]) - aGradR{2}/BetaVals(2)./circshift(mySqrtL, [0 -1 0]) + ...
                    aGradL{1}/BetaVals(1)./circshift(mySqrtR, [1 0 0]) + aGradL{2}/BetaVals(2)./circshift(mySqrtR, [0 1 0]);
    end
elseif ndims(toRegularize) == 3
    aGradL{1}=(toRegularize - circshift(toRegularize,[1 0 0]))*(epsR/BetaVals(1));  % cyclic rotation
    aGradL{2}=(toRegularize - circshift(toRegularize,[0 1 0]))*(epsR/BetaVals(2));  % cyclic rotation
    aGradL{3}=(toRegularize - circshift(toRegularize,[0 0 1]))*(epsR/BetaVals(3));  % cyclic rotation
    aGradR{1}=(circshift(toRegularize,[-1 0 0]) - toRegularize)*(epsR/BetaVals(1));
    aGradR{2}=(circshift(toRegularize,[0 -1 0]) - toRegularize)*(epsR/BetaVals(2));
    aGradR{3}=(circshift(toRegularize,[0 0 -1]) - toRegularize)*(epsR/BetaVals(3));
    mySqrtL = sqrt((abssqr(aGradL{1})+abssqr(aGradL{2})+abssqr(aGradL{3}))+(1.0-epsR)*toRegularizeSqr); 
    mySqrtR = sqrt((abssqr(aGradR{1})+abssqr(aGradR{2})+abssqr(aGradR{3}))+(1.0-epsR)*toRegularizeSqr);
    mySqrt = mySqrtL + mySqrtR; 
    mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
    myRegGrad = ((aGradL{1}*(epsR/BetaVals(1)) + aGradL{2}*(epsR/BetaVals(2)) + aGradL{3}*(epsR/BetaVals(3)))+ (1.0-epsR)*toRegularize)./mySqrtL - ...
                ((aGradR{1}*(epsR/BetaVals(1)) + aGradR{2}*(epsR/BetaVals(2)) + aGradR{3}*(epsR/BetaVals(3)))+ (1.0-epsR)*toRegularize)./mySqrtR - ...
                (aGradR{1}+ (1.0-epsR)*toRegularize)*(epsR/BetaVals(1))./circshift(mySqrtL, [-1 0 0]) - (aGradR{2}+ (1.0-epsR)*toRegularize)*(epsR/BetaVals(2))./circshift(mySqrtL, [0 -1 0]) - ...
                (aGradR{3}+ (1.0-epsR)*toRegularize)*(epsR/BetaVals(3))./circshift(mySqrtL, [0 0 -1]) + (aGradL{1}+ (1.0-epsR)*toRegularize)*(epsR/BetaVals(1))./circshift(mySqrtR, [1 0 0]) + ...
                (aGradL{2}+ (1.0-epsR)*toRegularize)*(epsR/BetaVals(2))./circshift(mySqrtR, [0 1 0]) + (aGradL{3}+ (1.0-epsR)*toRegularize)*(epsR/BetaVals(3))./circshift(mySqrtR, [0 0 1]) ;
else % 1-D
    aGradL=(toRegularize - circshift(toRegularize,1))/BetaVals(1);  
    aGradR=(circshift(toRegularize,-1) - toRegularize)/BetaVals(1);  % cyclic rotation
    mySqrtL = sqrt(epsR*abssqr(aGradL)+(1.0-epsR)*toRegularizeSqr);
    mySqrtR = sqrt(epsR*abssqr(aGradR)+(1.0-epsR)*toRegularizeSqr);
    mySqrt = mySqrtL + mySqrtR; 
    mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
    myRegGrad = (epsR*aGradL/BetaVals(1)+ (1.0-epsR)*2*toRegularize)./mySqrtL - (epsR*aGradR/BetaVals(1)+ (1.0-epsR)*2*toRegularize)./mySqrtR - ...
                aGradR/BetaVals(1)./circshift(mySqrtL, -1) + aGradL/BetaVals(1)./circshift(mySqrtR, 1);
end
myReg = sum(mySqrt);

%if ~isreal(myRegGrad)
%    myRegGrad=conj(myRegGrad);
%end