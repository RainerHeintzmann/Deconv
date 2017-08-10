function myGrad=BwdHalfResidIlluObjConvPSF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
myrft=residuum .* conj(myOtf);  % residuum is already in Fourier space due to the Ptychography trick!
global aResampling;
if any(aResampling~=1)
       myrft=rft_resize(myrft,aResampling);  % if the user wants to use a different reconstruction grid
       norm3D=norm3D/ prod(aResampling);
end

myGrad = norm3D * aRecon .* rift(myrft);  