function Recons=FwdObjIlluConvPSF(aRecon,ftRecons,myIllum,myOtf,norm3D)
ftRecons=rft(aRecon .* myIllum);
Recons=norm3D*rift(ftRecons .* myOtf);  % convolve with corresponding PSF, still in Fourier space

