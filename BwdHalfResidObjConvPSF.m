function myGrad=BwdHalfResidObjConvPSF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D) % assumes that the residuum already lives in Fourier space. This saves time for Gaussian noise models
% It is also what was labeled "Fourier-Ptychography"
myrft=residuum .* conj(myOtf);
global aResampling;
if any(aResampling~=1)
       myrft=rft_resize(myrft,aResampling);  % if the user wants to use a different reconstruction grid
       norm3D=norm3D/ prod(aResampling);
end

myGrad = norm3D*rift(myrft);
