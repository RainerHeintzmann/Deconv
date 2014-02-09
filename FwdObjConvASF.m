function Recons=FwdObjConvASF(aRecon,ftRecons,myIllum,myOtf,norm3D)
Recons=norm3D*ift(ftRecons .* myOtf);  % convolve with corresponding PSF, still in Fourier space

