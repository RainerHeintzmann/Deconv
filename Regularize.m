function [myReg,myRegGrad]=Regularize(toRegularize,BetaVals)
global RegularisationParameters;

myReg=0;
myRegGrad=0;
 
myLambda=RegularisationParameters(1,1); %case 'GS'  % Gradient roughness penalty: Gradient ^ 2
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGS(toRegularize,BetaVals);
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(2,1); %case 'AR'  % Arigovindan's roughness penalty
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeAR(toRegularize,BetaVals,RegularisationParameters(2,2));
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
    [aReg,aRegGrad]=RegularizeGR(toRegularize,BetaVals,RegularisationParameters(5,2),RegularisationParameters(5,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end

myLambda=RegularisationParameters(11,1); %case 'CO'  % Conchello's penalty based on the absolute sqr of the object
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeCO(toRegularize);  % Conchello regularisation g^2
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
