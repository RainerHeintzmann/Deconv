function Recons=FwdObjIlluConvPSF(aRecon,ftRecons,myIllum,myOtf,norm3D)
global aResampling
ftRecons=rft(aRecon .* myIllum);
if any(aResampling~=1)
       ftRecons=rft_resize(ftRecons,1./aResampling);  % if the user wants to use a different reconstruction grid
end

Recons=norm3D*rift(ftRecons .* myOtf);  % convolve with corresponding PSF, still in Fourier space

