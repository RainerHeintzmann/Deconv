function myGrad=BwdResidFT(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
global ComplexObj;
if ~ComplexObj
    myGrad = real(norm3D*ift(residuum));
else
    myGrad = norm3D*ift(residuum);
end