% [myReg,myRegGrad]=RegularizeTV(toRegularize,BetaVals,epsR) computes Total variation regularisation
% penalty = sqrt(|grad(f)|^2+epsR)
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient
% The total variation code below is based on a simple 2-point finite difference
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

function [myReg,myRegGrad]=RegularizeTV(toRegularize,BetaVals,epsR)
epsC=1e-10;
if nargin < 3
    epsR=1e2*1e2;  % Modifies the TV norm to behave better (e.g. considering digitalisation effects)
end
if (ndims(toRegularize) == 2) || (size(toRegularize,3) == 1)
    if (ndims(toRegularize) == 2)
        aGradL{1}=(toRegularize - circshift(toRegularize,[1 0]))/BetaVals(1);  % cyclic rotation
        aGradL{2}=(toRegularize - circshift(toRegularize,[0 1]))/BetaVals(2);  % cyclic rotation
        aGradR{1}=(circshift(toRegularize,[-1 0]) - toRegularize)/BetaVals(1);
        aGradR{2}=(circshift(toRegularize,[0 -1]) - toRegularize)/BetaVals(2);
        mySqrtL = sqrt(abssqr(aGradL{1}) + abssqr(aGradL{2}) + epsR);
        mySqrtR = sqrt(abssqr(aGradR{1}) + abssqr(aGradR{2}) + epsR);
        mySqrt = mySqrtL + mySqrtR;
        mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
        % gradient is tested
        myRegGrad = (aGradL{1}/BetaVals(1) + aGradL{2}/BetaVals(2))./mySqrtL - (aGradR{1}/BetaVals(1) + aGradR{2}/BetaVals(2))./mySqrtR - ...
                    aGradR{1}/BetaVals(1)./circshift(mySqrtL, [-1 0]) - aGradR{2}/BetaVals(2)./circshift(mySqrtL, [0 -1]) + ...
                    aGradL{1}/BetaVals(1)./circshift(mySqrtR, [1 0]) + aGradL{2}/BetaVals(2)./circshift(mySqrtR, [0 1]);
    else (size(toRegularize,3) == 1)
        aGradL{1}=(toRegularize - circshift(toRegularize,[1 0 0]))/BetaVals(1);  % cyclic rotation
        aGradL{2}=(toRegularize - circshift(toRegularize,[0 1 0]))/BetaVals(2);  % cyclic rotation
        aGradR{1}=(circshift(toRegularize,[-1 0 0]) - toRegularize)/BetaVals(1);
        aGradR{2}=(circshift(toRegularize,[0 -1 0]) - toRegularize)/BetaVals(2);
        mySqrtL = sqrt(abssqr(aGradL{1})+abssqr(aGradL{2})+epsR);
        mySqrtR = sqrt(abssqr(aGradR{1})+abssqr(aGradR{2})+epsR);
        mySqrt =  mySqrtL + mySqrtR;
        mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
        myRegGrad = (aGradL{1}/BetaVals(1) + aGradL{2}/BetaVals(2))./mySqrtL - (aGradR{1}/BetaVals(1) + aGradR{2}/BetaVals(2))./mySqrtR - ...
                    aGradR{1}/BetaVals(1)./circshift(mySqrtL, [-1 0 0]) - aGradR{2}/BetaVals(2)./circshift(mySqrtL, [0 -1 0]) + ...
                    aGradL{1}/BetaVals(1)./circshift(mySqrtR, [1 0 0]) + aGradL{2}/BetaVals(2)./circshift(mySqrtR, [0 1 0]);
    end
elseif ndims(toRegularize) == 3
    aGradL{1}=(toRegularize - circshift(toRegularize,[1 0 0]))/BetaVals(1);  % cyclic rotation
    aGradL{2}=(toRegularize - circshift(toRegularize,[0 1 0]))/BetaVals(2);  % cyclic rotation
    aGradL{3}=(toRegularize - circshift(toRegularize,[0 0 1]))/BetaVals(3);  % cyclic rotation
    aGradR{1}=(circshift(toRegularize,[-1 0 0]) - toRegularize)/BetaVals(1);
    aGradR{2}=(circshift(toRegularize,[0 -1 0]) - toRegularize)/BetaVals(2);
    aGradR{3}=(circshift(toRegularize,[0 0 -1]) - toRegularize)/BetaVals(3);
    mySqrtL = sqrt(abssqr(aGradL{1})+abssqr(aGradL{2})+abssqr(aGradL{3})+epsR); 
    mySqrtR = sqrt(abssqr(aGradR{1})+abssqr(aGradR{2})+abssqr(aGradR{3})+epsR);
    mySqrt = mySqrtL + mySqrtR; 
    mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
    myRegGrad = (aGradL{1}/BetaVals(1) + aGradL{2}/BetaVals(2) + aGradL{3}/BetaVals(3))./mySqrtL - ...
                (aGradR{1}/BetaVals(1) + aGradR{2}/BetaVals(2) + aGradR{3}/BetaVals(3))./mySqrtR - ...
                aGradR{1}/BetaVals(1)./circshift(mySqrtL, [-1 0 0]) - aGradR{2}/BetaVals(2)./circshift(mySqrtL, [0 -1 0]) - ...
                aGradR{3}/BetaVals(3)./circshift(mySqrtL, [0 0 -1]) + aGradL{1}/BetaVals(1)./circshift(mySqrtR, [1 0 0]) + ...
                aGradL{2}/BetaVals(2)./circshift(mySqrtR, [0 1 0]) + aGradL{3}/BetaVals(3)./circshift(mySqrtR, [0 0 1]);
else % 1-D
    aGradL=(toRegularize - circshift(toRegularize,1))/BetaVals(1);  
    aGradR=(circshift(toRegularize,-1) - toRegularize)/BetaVals(1);  % cyclic rotation
    mySqrtL = sqrt(abssqr(aGradL)+epsR);
    mySqrtR = sqrt(abssqr(aGradR)+epsR);
    mySqrt = mySqrtL + mySqrtR; 
    mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
    myRegGrad = aGradL/BetaVals(1)./mySqrtL - aGradR/BetaVals(1)./mySqrtR - ...
                aGradR/BetaVals(1)./circshift(mySqrtL, -1) + aGradL/BetaVals(1)./circshift(mySqrtR, 1);
end
myReg = sum(mySqrt);

%if ~isreal(myRegGrad)
%    myRegGrad=conj(myRegGrad);
%end