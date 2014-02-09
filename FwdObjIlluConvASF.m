function Recons=FwdObjIlluConvASF(aRecon,ftRecons,myIllum,myOtf,norm3D)
ftRecons=ft(aRecon .* myIllum);
Recons=norm3D*ift(ftRecons .* myOtf);  % convolve with corresponding PSF, still in Fourier space

