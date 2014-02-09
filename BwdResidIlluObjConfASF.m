function myGrad=BwdResidIlluObjConfASF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
global ComplexObj;
myGrad = aRecon .* norm3D*ift(ft(residuum) .* conj(myOtf));

if ~ComplexObj
    myGrad = real(aRecon .* norm3D*ift(ft(residuum) .* conj(myOtf)));
else
    myGrad = aRecon .* norm3D*ift(ft(residuum) .* conj(myOtf));
end
