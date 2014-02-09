function res=ApplyBgModel(aRecon2,ftRecons,myIllum,myOtf,norm3D,afkt,BgValue)
res=afkt(aRecon2,ftRecons,myIllum,myOtf,norm3D)+BgValue; % The Bg Model is unaffected

