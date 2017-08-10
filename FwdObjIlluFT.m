function Recons=FwdObjIlluFT(aRecon,ftRecons,myIllum,myOtf,norm3D)
Recons=ft(aRecon .* myIllum);
%Recons=norm3D*ftRecons;  % convolve with corresponding PSF, still in Fourier space

