function myGrad=BwdResidIlluObjConfPSF_Slice(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
myrft=repmat(rft(residuum),[1 1 size(myOtf,3)]) .* conj(myOtf);

global aResampling;
if any(aResampling~=1)
       myrft=rft_resize(myrft,aResampling);  % if the user wants to use a different reconstruction grid
end

myGrad = norm3D * aRecon .* rift(myrft);  