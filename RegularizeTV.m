% [myReg,myRegGrad]=RegularizeTV(toRegularize,BetaVals,epsR) computes Total variation regularisation
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient
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

function [myReg,myRegGrad]=RegularizeTV(toRegularize,BetaVals,epsR)
epsC=1e-10;
if nargin < 3
    epsR=1e2*1e2;  % Modifies the TV norm to behave better (e.g. considering digitalisation effects)
end
if (ndims(toRegularize) == 2) || (size(toRegularize,3) == 1)
    if (ndims(toRegularize) == 2)
        aGrad{1}=(circshift(toRegularize,[1 0])-circshift(toRegularize,[-1 0]))/(2*BetaVals(1));  % cyclic rotation
        aGrad{2}=(circshift(toRegularize,[0 1])-circshift(toRegularize,[0 -1]))/(2*BetaVals(2));  % cyclic rotation
        mySqrt = sqrt(abssqr(aGrad{1})+abssqr(aGrad{2})+epsR);
        mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
        myRegGrad = (1/(4*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,[-2 0]))./circshift(mySqrt,[-1 0]) - (circshift(toRegularize,[2 0])-toRegularize)./circshift(mySqrt,[1 0])) + ...
            (1/(4*BetaVals(2)*BetaVals(2)))*((toRegularize-circshift(toRegularize,[0 -2]))./circshift(mySqrt,[0 -1]) - (circshift(toRegularize,[0 2])-toRegularize)./circshift(mySqrt,[0 1]));
    else (size(toRegularize,3) == 1)
        aGrad{1}=(circshift(toRegularize,[1 0 0])-circshift(toRegularize,[-1 0 0]))/(2*BetaVals(1));  % cyclic rotation
        aGrad{2}=(circshift(toRegularize,[0 1 0])-circshift(toRegularize,[0 -1 0]))/(2*BetaVals(2));  % cyclic rotation
        mySqrt = sqrt(abssqr(aGrad{1})+abssqr(aGrad{2})+epsR);
        mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
        myRegGrad = (1/(4*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,[-2 0 0]))./circshift(mySqrt,[-1 0 0]) - (circshift(toRegularize,[2 0 0])-toRegularize)./circshift(mySqrt,[1 0 0])) + ...
            (1/(4*BetaVals(2)*BetaVals(2)))*((toRegularize-circshift(toRegularize,[0 -2 0]))./circshift(mySqrt,[0 -1 0]) - (circshift(toRegularize,[0 2 0])-toRegularize)./circshift(mySqrt,[0 1 0]));
    end
elseif ndims(toRegularize) == 3
    aGrad{1}=(circshift(toRegularize,[1 0 0])-circshift(toRegularize,[-1 0 0]))/(2*BetaVals(1));  % cyclic rotation
    aGrad{2}=(circshift(toRegularize,[0 1 0])-circshift(toRegularize,[0 -1 0]))/(2*BetaVals(2));  % cyclic rotation
    aGrad{3}=(circshift(toRegularize,[0 0 1])-circshift(toRegularize,[0 0 -1]))/(2*BetaVals(3));  % cyclic rotation
    mySqrt = sqrt(abssqr(aGrad{1})+abssqr(aGrad{2})+abssqr(aGrad{3})+epsR);
    mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
    myRegGrad = (1/(4*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,[-2 0 0]))./circshift(mySqrt,[-1 0 0]) - (circshift(toRegularize,[2 0 0])-toRegularize)./circshift(mySqrt,[1 0 0])) + ...
        (1/(4*BetaVals(2)*BetaVals(2)))*((toRegularize-circshift(toRegularize,[0 -2 0]))./circshift(mySqrt,[0 -1 0]) - (circshift(toRegularize,[0 2 0])-toRegularize)./circshift(mySqrt,[0 1 0])) + ...
        (1/(4*BetaVals(3)*BetaVals(3)))*((toRegularize-circshift(toRegularize,[0 0 -2]))./circshift(mySqrt,[0 0 -1]) - (circshift(toRegularize,[0 0 2])-toRegularize)./circshift(mySqrt,[0 0 1]));
else % 1-D
    aGrad{1}=(circshift(toRegularize,1)-circshift(toRegularize,-1))/(2*BetaVals(1));  % cyclic rotation
    mySqrt = sqrt(abssqr(aGrad{1})+epsR);
    mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
    myRegGrad = (1/(4*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,-2))./circshift(mySqrt,-1) - (circshift(toRegularize,2)-toRegularize)./circshift(mySqrt,1));
end
myReg = sum(mySqrt);

%if ~isreal(myRegGrad)
%    myRegGrad=conj(myRegGrad);
%end