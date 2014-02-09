function myGrad=BwdResidObjConfASFSqr(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)
%global ComplexObj;
global ReconsSaved;  % Was set in the forward model

%if ~ComplexObj
%    myGrad = real(norm3D*ift(ft(residuum) .* conj(myOtf)));
%else
    myGrad = norm3D*ift(ft(residuum) .* conj(myOtf));
%end

residuum = 2*ReconsSaved.*residuum;  % This goes back from the intensity world to the amplitude world
clear ReconsSaved;
