function Recons=FwdObjFT(aRecon,ftRecons,myIllum,myOtf,norm3D)
Recons=norm3D*ftRecons;  % convolve with corresponding PSF, still in Fourier space

