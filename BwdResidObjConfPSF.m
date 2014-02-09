function myGrad=BwdResidObjConfPSF(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
myrft=rft(residuum) .* conj(myOtf);
global aResampling;
if any(aResampling~=1)
       myrft=rft_resize(myrft,aResampling);  % if the user wants to use a different reconstruction grid
end

myGrad = norm3D*rift(myrft);
