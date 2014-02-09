function myGrad=BwdResidObjIlluConfPSF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
myGrad = myIllum .* norm3D * rift(rft(residuum) .* conj(myOtf));
