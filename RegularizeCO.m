% [myReg,myRegGrad]=RegularizeCO(toRegularize) computes Conchello's penalty based on object^2
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient
% 
% See Footnote in: Conchello et al. Science 263, 1483-1487, 1995 

function [myReg,myRegGrad]=RegularizeCO(toRegularize)
myReg = sum(abssqr(toRegularize)); 
myRegGrad = 2*toRegularize;
