function Recons=FwdObjIlluConvPSF_Slice(aRecon,ftRecons,myIllum,myOtf,norm3D)  % Performs the forward convolution but also extracting only the middle slice
global aResampling
ftRecons=rft(aRecon .* myIllum);
% OTF was modifiend in GenericDeconvolution the be consisten with summing in Fourier space for RFT:
% myOtf=myOtf .* ifftshift(exp(i*(floor(size(myOtf,3)/2)/size(myOtf,3))*2*pi*zz(size(myOtf)))) /sqrt(size(myOtf,3));  % modifies the OTF the be consisten with summing in Fourier space for RFT
if any(aResampling~=1)
       ftRecons=rft_resize(ftRecons,1./aResampling);  % if the user wants to use a different reconstruction grid
end

Recons=norm3D*rift(sum(ftRecons .* myOtf,[],3));  % convolve with corresponding PSF, still in Fourier space but extract only the central slice

