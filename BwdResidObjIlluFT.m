function myGrad=BwdResidObjIlluFT(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
global ComplexObj;
if ~ComplexObj
    myGrad = real(aRecon .* ift(residuum)); % norm3D*
else
    myGrad = aRecon .* ift(residuum); % norm3D*
end
