function Recons=FwdObjConvPSF(aRecon2,ftRecons,myIllum,myOtf,norm3D)

Recons=norm3D*rift(ftRecons .* myOtf);  % convolve with corresponding PSF, still in Fourier space

global aRecon;
clear aRecon;  % is not needed any longer after the forward convolution
