function Recons = FwdNormMeasSum(aRecon2,ftRecons,myIllum,myOtf,norm3D,FwdModel)
global measSum;
global FwdSum;
global NormRecons;
global DeconvMask;
Recons = FwdModel(aRecon2,ftRecons,myIllum,myOtf,norm3D);
FwdSum = sum(Recons,DeconvMask); 
Recons = Recons*measSum/FwdSum;
NormRecons=Recons;
