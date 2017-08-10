function myGrad=BwdResidObjConvASF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
global ComplexObj;
if ~ComplexObj
    myGrad = real(norm3D*ift(ft(residuum) .* conj(myOtf)));
else
    myGrad = norm3D*ift(ft(residuum) .* conj(myOtf));
end