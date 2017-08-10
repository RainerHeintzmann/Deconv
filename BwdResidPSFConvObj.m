function myGrad=BwdResidPSFConvObj(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)  % This is used for blind PSF deconvolutions
%myGrad=norm3D*rift(rft(residuum/sum(aRecon)) .* conj(ftRecon));
myrft=rft(residuum);
global aResampling;
if any(aResampling~=1)
       myrft=rft_resize(myrft,aResampling);  % if the user wants to use a different reconstruction grid
       norm3D=norm3D/ prod(aResampling);
end

myGrad=norm3D*rift(myrft .* conj(ftRecon));
