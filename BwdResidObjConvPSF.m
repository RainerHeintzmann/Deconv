function myGrad=BwdResidObjConvPSF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
global aResampling;
myrft=rft(residuum) .* conj(myOtf);
if any(aResampling~=1)
       myrft=rft_resize(myrft,aResampling);  % if the user wants to use a different reconstruction grid
       norm3D=norm3D/ prod(aResampling);
end

myGrad = norm3D*rift(myrft);
