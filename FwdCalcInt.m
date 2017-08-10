function Recons = FwdCalcInt(aRecon,ftRecon,myIllum,myOtf,norm3D, myfun)
global ampRecon;
Recons = myfun(aRecon,ftRecon,myIllum,myOtf,norm3D);
ampRecon = Recons;

Recons = abssqr(Recons);