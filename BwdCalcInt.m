function myGrad = BwdCalcInt(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D,myfun)
global ampRecon;
myGrad = 2*real(residuum).*ampRecon;

myGrad = myfun(myGrad,aRecon,ftRecon,myIllum,myOtf,norm3D);