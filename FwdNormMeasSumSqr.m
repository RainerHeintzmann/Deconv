function Recons = FwdNormMeasSumSqr(aRecon2,ftRecons,myIllum,myOtf,norm3D,FwdModel)
global measSumSqr;
global FwdSumSqr;
global NormRecons;
global DeconvMask;
Recons = FwdModel(aRecon2,ftRecons,myIllum,myOtf,norm3D);
% FwdSumSqr = sum(abs(Recons).^2);
FwdSumSqr = sqrt(sum(abs(Recons).^2,DeconvMask));
Recons = Recons*measSumSqr/FwdSumSqr;
NormRecons=Recons;
