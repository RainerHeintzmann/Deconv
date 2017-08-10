function myGrad=BwdResidObjIlluConvASF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
global ComplexObj;
if ~ComplexObj
    if isreal(myIllum)
        myGrad = real(myIllum .* norm3D*ift(ft(residuum) .* conj(myOtf)));
    else
        myGrad = real(conj(myIllum) .* norm3D*ift(ft(residuum) .* conj(myOtf)));
    end
else
    if isreal(myIllum)
        myGrad = myIllum .* norm3D*ift(ft(residuum) .* conj(myOtf));
    else
        myGrad = conj(myIllum) .* norm3D*ift(ft(residuum) .* conj(myOtf));
    end
end
