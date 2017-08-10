% [myReg,myRegGrad]=RegularizeNegSqr(toRegularize) computes Negative Value penalty regularisation
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient

function [myReg,myRegGrad]=RegularizeNegSqr(toRegularize)
myReg = sum(abssqr(toRegularize).*(toRegularize<0));  % Just affects the real part
myRegGrad = 2*toRegularize .* (toRegularize<0);
