function myGrad = BwdNormMeasSumSqr(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D,BwdModel)
global measSumSqr;
global FwdSumSqr;
global NormRecons;

% residuum = convertCpxObjToVec(residuum);
% normImg = convertCpxObjToVec(NormRecons);
        
%     residuum1 = residuum(1:end/2);
%     residuum2 = residuum(end/2+1:end);
%     normImg1 = normImg(1:end/2);
%     normImg2 = normImg(end/2+1:end);
residuumR = real(residuum);
residuumI = imag(residuum);
NormReconsR = real(NormRecons);
NormReconsI = imag(NormRecons);

%     derivRR = measSumSqr*residuum1'/FwdSumSqr - 2*(residuum1'*normImg1)*normImg1'/measSumSqr;
%     derivIR = -2*residuum2'*normImg2*normImg1'/measSumSqr;
%     derivRI = -2*residuum1'*normImg1*normImg2'/measSumSqr;
%     derivII = measSumSqr*residuum2'/FwdSumSqr - 2*(residuum2'*normImg2)*normImg2'/measSumSqr;

% derivRR = measSumSqr*residuumR/FwdSumSqr - 2*sum(residuumR.*NormReconsR)*NormReconsR/measSumSqr;
% derivIR = -2*sum(residuumI.*NormReconsI)*NormReconsR/measSumSqr;
% derivRI = -2*sum(residuumR.*NormReconsR)*NormReconsI/measSumSqr;
% derivII = measSumSqr*residuumI/FwdSumSqr - 2*sum(residuumI.*NormReconsI)*NormReconsI/measSumSqr;

derivRR = measSumSqr*residuumR - sum(residuumR.*NormReconsR)*NormReconsR/measSumSqr;
derivIR = -sum(residuumI.*NormReconsI)*NormReconsR/measSumSqr;
derivRI = -sum(residuumR.*NormReconsR)*NormReconsI/measSumSqr;
derivII = measSumSqr*residuumI - sum(residuumI.*NormReconsI)*NormReconsI/measSumSqr;

% myGrad = [derivRR + derivIR, derivRI + derivII]/FwdSumSqr;
% myGrad = derivRR + derivIR + 1i*(derivRI + derivII);
myGrad = (derivRR + derivIR + 1i*(derivRI + derivII))/FwdSumSqr;
   
% myGrad = convertVecToCpxObj(myGrad);
    
myGrad = BwdModel(myGrad,aRecon,ftRecon,myIllum,myOtf,norm3D);
