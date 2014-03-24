function myGrad = BwdNormMeasSum(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D,BwdModel)
global measSum;
global FwdSum;
global NormRecons;

% residuum = convertObjToVec(residuum);
% normImg = convertObjToVec(NormRecons);
    
myGrad = (measSum*residuum - sum(residuum.*NormRecons))/FwdSum; % (residuum'*measSum-residuum'*normImg)/FwdSum;
% myGrad = convertVecToObj(myGrad);
    
myGrad = BwdModel(myGrad,aRecon,ftRecon,myIllum,myOtf,norm3D);
