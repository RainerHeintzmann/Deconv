function Recons=FwdHalfObjConvPSF(aRecon2,ftRecons,myIllum,myOtf,norm3D)  % will lead to a Fwd projection that still lives in Fourier space. This saves time for Gaussian noise models
% It is also what was labeled "Fourier-Ptychography"

Recons=norm3D * ftRecons .* myOtf;  % convolve with corresponding PSF, but stay in Fourier space

global aRecon;
clear aRecon;  % is not needed any longer after the forward convolution
