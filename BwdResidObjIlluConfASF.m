function myGrad=BwdResidObjIlluConfASF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
global ComplexObj;
if ~ComplexObj
    myGrad = real(myIllum .* norm3D*ift(ft(residuum) .* conj(myOtf)));
else
    myGrad = myIllum .* norm3D*ift(ft(residuum) .* conj(myOtf));
end
