function mybool=NeedsRegularisation()
global RegularisationParameters;
AreReg=[1 2 3 4 5 11];  % These parameter positions correspond to lambdas in the regularisation
mybool=any(RegularisationParameters(AreReg,1));

