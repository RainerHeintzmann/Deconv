function myGrad=BwdResidIlluObjConfPSF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
myGrad = aRecon .* norm3D * rift(rft(residuum) .* conj(myOtf));
