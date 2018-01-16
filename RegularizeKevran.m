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
    aGrad{1}=(circshift(toRegularize,[-1 0 0]) - circshift(toRegularize,[1 0 0]))*(epsR/BetaVals(1)/2);  % cyclic rotation
    aGrad{2}=(circshift(toRegularize,[0 -1 0]) - circshift(toRegularize,[0 1 0]))*(epsR/BetaVals(2)/2);  % cyclic rotation
    mySqrt = sqrt(abssqr(aGrad{1})+abssqr(aGrad{2})+(1.0-epsR)*toRegularizeSqr); 
    mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
    tmp = toRegularize./mySqrt;
    myRegGrad = (1-epsR)*tmp + (epsR/BetaVals(1)/2)*(circshift(aGrad{1} ./ mySqrt,[1 0 0]) - circshift(aGrad{1} ./ mySqrt,[-1 0 0])) + ...
    (epsR/BetaVals(2)/2).*(circshift(aGrad{2} ./ mySqrt,[0 1 0]) - circshift(aGrad{2} ./ mySqrt,[0 -1 0]));

elseif ndims(toRegularize) == 3
    aGrad{1}=(circshift(toRegularize,[-1 0 0]) - circshift(toRegularize,[1 0 0]))*(epsR/BetaVals(1)/2);  % cyclic rotation
    aGrad{2}=(circshift(toRegularize,[0 -1 0]) - circshift(toRegularize,[0 1 0]))*(epsR/BetaVals(2)/2);  % cyclic rotation
    aGrad{3}=(circshift(toRegularize,[0 0 -1]) - circshift(toRegularize,[0 0 1]))*(epsR/BetaVals(3)/2);  % cyclic rotation
    mySqrt = sqrt(abssqr(aGrad{1})+abssqr(aGrad{2})+abssqr(aGrad{3})+(1.0-epsR)*toRegularizeSqr); 
    mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
    tmp = toRegularize./mySqrt;
    myRegGrad = (1-epsR)*tmp + (epsR/BetaVals(1)/2)*(circshift(aGrad{1} ./ mySqrt,[1 0 0]) - circshift(aGrad{1} ./ mySqrt,[-1 0 0])) + ...
    (epsR/BetaVals(2)/2).*(circshift(aGrad{2} ./ mySqrt,[0 1 0]) - circshift(aGrad{2} ./ mySqrt,[0 -1 0])) + ...
    (epsR/BetaVals(3)/2).*(circshift(aGrad{3} ./ mySqrt,[0 0 1]) - circshift(aGrad{3} ./ mySqrt,[0 0 -1])) ;
else % 1-D
    aGrad{1}=(circshift(toRegularize,[-1 0 0]) - circshift(toRegularize,[1 0 0]))*(epsR/BetaVals(1)/2);  % cyclic rotation
    mySqrt = sqrt(abssqr(aGrad{1})+(1.0-epsR)*toRegularizeSqr); 
    mySqrt(mySqrt< epsC) = epsC;  % To avoid divisions by zero
    tmp = toRegularize./mySqrt;
    myRegGrad = (1-epsR)*tmp + (epsR/BetaVals(1)/2)*(circshift(aGrad{1} ./ mySqrt,[1 0 0]) - circshift(aGrad{1} ./ mySqrt,[-1 0 0]));
end
myReg = sum(mySqrt);

%if ~isreal(myRegGrad)
%    myRegGrad=conj(myRegGrad);
%end