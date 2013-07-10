% [myReg,myRegGrad]=RegularizeNegSqr(toRegularize) computes Negative Value penalty regularisation
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient

function [myReg,myRegGrad]=RegularizeTV(toRegularize)
myReg = sum(toRegularize.^2.*(toRegularize<0));
myRegGrad = 2*toRegularize .* (toRegularize<0);
