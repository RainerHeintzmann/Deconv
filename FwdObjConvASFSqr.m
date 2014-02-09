function Recons=FwdObjConvASFSqr(aRecon,ftRecons,myIllum,myOtf,norm3D)
global ReconsSaved;  % This has to be remembered, since it will be needed for the backward model

Recons=norm3D*ift(ftRecons .* myOtf);  % convolve with corresponding PSF, still in Fourier space

ReconsSaved=Recons;
Recons=abssqr(Recons);  % After convolution the intensity is now calculated to be compared to the intensity data
