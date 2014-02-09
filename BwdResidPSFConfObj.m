function myGrad=BwdResidPSFConfObj(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)  % This is used for blind PSF deconvolutions
myGrad = norm3D*rift(rft(fftshift(residuum/sum(aRecon))) .* conj(ftRecon));
