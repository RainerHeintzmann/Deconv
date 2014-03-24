function myGrad=BwdResidPSFConfObj(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)  % This is used for blind PSF deconvolutions
%myGrad=norm3D*rift(rft(residuum/sum(aRecon)) .* conj(ftRecon));
myGrad=norm3D*rift(rft(residuum) .* conj(ftRecon));
