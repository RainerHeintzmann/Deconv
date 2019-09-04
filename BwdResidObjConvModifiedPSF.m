function myGrad=BwdResidObjConvModifiedPSF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
global BwdOTF;
global aResampling;
myrft=rft(residuum) .* BwdOTF;
if any(aResampling~=1)
       myrft=rft_resize(myrft,aResampling);  % if the user wants to use a different reconstruction grid
       norm3D=norm3D/ prod(aResampling);
end

myGrad = norm3D*rift(myrft);
