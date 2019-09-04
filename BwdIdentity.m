function myGrad=BwdIdentity(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
global aResampling;
if any(aResampling~=1)
    myrft=rft(residuum);
    myrft=rft_resize(myrft,aResampling);  % if the user wants to use a different reconstruction grid
    norm3D=norm3D/ prod(aResampling);
    myGrad = norm3D*rift(myrft);
else
    myGrad = residuum;
end

