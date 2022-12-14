function [myReg,myRegGrad]=Regularize(toRegularize,BetaVals)
global RegularisationParameters;
global RefImgX;
global RefImgY;

myReg=0;
myRegGrad=0;
 
myLambda=RegularisationParameters(1,1); %case 'GS'  % Gradient roughness penalty: Gradient ^ 2
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGS(toRegularize,BetaVals);
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(2,1); %case 'AR'  % Arigovindan's roughness penalty
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeER(toRegularize,BetaVals,RegularisationParameters(2,2));
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

myLambda=RegularisationParameters(22,1); %case 'LAP' Laplaceoperator 6 star scheme  
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeLap(toRegularize,BetaVals); 
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(23,1); %case 'GRLapGrad6'  % Good's roughness Gradient penalty: Laplaceoperator ^ 2/f (central difference 6 star scheme)
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGRLapGrad6(toRegularize,BetaVals,RegularisationParameters(23,2),RegularisationParameters(23,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(24,1); %case 'GRLapGradReg'  % Good's roughness Gradient + current value penalty: Laplaceoperator ^ 2/f (6 star scheme)
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGRLapGradReg(toRegularize,BetaVals,RegularisationParameters(24,2),RegularisationParameters(24,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(25,1); %case 'GRCentral'  % Good's roughness penalty: Gradient ^ 2/f (central differences)
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGRZentrale(toRegularize,BetaVals,RegularisationParameters(25,2),RegularisationParameters(25,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(26,1); %case 'GRLap6'  % Good's roughness penalty: Laplaceopertor /f
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGRLap6(toRegularize,BetaVals,RegularisationParameters(26,2),RegularisationParameters(26,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(27,1); %case 'Lap27'  % Laplaceoperator 27star scheme
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeLap27(toRegularize,BetaVals,RegularisationParameters(27,2),RegularisationParameters(27,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(28,1); %case 'GRLapGradReg'  % Good's roughness Gradient + current value penalty: Laplaceoperator ^ 2/f (6 star scheme)
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGRLapDivReg(toRegularize,BetaVals,RegularisationParameters(28,2),RegularisationParameters(28,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end
myLambda=RegularisationParameters(29,1); %case 'GRLapGradReg'  % Good's roughness Gradient + current value penalty: Laplaceoperator ^ 2/f (6 star scheme)
if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGRLapDiv6(toRegularize,BetaVals,RegularisationParameters(29,2),RegularisationParameters(29,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end


myLambda=RegularisationParameters(30,1); %case 'GRLapGradReg'  % Good's roughness Gradient + current value penalty: Laplaceoperator ^ 2/f (6 star scheme)

if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGRLapDiv6stabil(toRegularize,BetaVals,RegularisationParameters(30,2),RegularisationParameters(30,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end

myLambda=RegularisationParameters(31,1); %case 'GRLapGradReg'  % Good's roughness Gradient + current value penalty: Laplaceoperator ^ 2/f (6 star scheme)

if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeGRLapDivReg2(toRegularize,BetaVals,RegularisationParameters(31,2),RegularisationParameters(31,3));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end

myLambda=RegularisationParameters(32,1); %case 'Kevran'  % Good's roughness Gradient + current value penalty: Laplaceoperator ^ 2/f (6 star scheme)

if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeKevran(toRegularize,BetaVals,RegularisationParameters(32,2));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end


myLambda=RegularisationParameters(33,1); %case CLE-GS with Guidance image

if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeCLE_GS(toRegularize,BetaVals,RefImgX,RefImgY,RegularisationParameters(33,2));%fengjiao 25.05.2020
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end

myLambda=RegularisationParameters(36,1); %case IG with Guidance image

if myLambda ~= 0
    [aReg,aRegGrad]=RegularizeIG(toRegularize,BetaVals,RefImgX,RegularisationParameters(36,2));
    myReg = myReg+myLambda * aReg; myRegGrad = myRegGrad + myLambda * aRegGrad;
end

end