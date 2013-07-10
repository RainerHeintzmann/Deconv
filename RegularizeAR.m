% [myReg,myRegGrad]=RegularizeAR(toRegularize,BetaVals) computes Arigovindan's roughness regularisation
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient
% This function is based on a talk by 
% Muthuvel  Arigovindan*,  Daniel Elnatan,  Jennifer  Fung,  Eric Branlund,  John W. Sedat,   and  David A.  Agard
% University of California at San Francisco
% suggesting a regolarisation term as given below

function [myReg,myRegGrad]=RegularizeAR(toRegularize,BetaVals)
        alpha= 1.0;  % Weighting between laplacian and noise penalty
        %alpha= 10000.0;  % Weighting between laplacian and noise penalty
        gamma = 1.0; % Weight for the intensity penalty under the logarithm
        aReconSqr=toRegularize.*toRegularize;

        if ndims(toRegularize) == 2 || (size(toRegularize,3) == 1)
            s1=1/1.6;s2=0.1727;
            H11=dip_convolve1d(toRegularize,s1*[1 -2 1],0,1)/(BetaVals(1)*BetaVals(1)); % second X derivative
            H22=dip_convolve1d(toRegularize,s1*[1 -2 1],1,1)/(BetaVals(2)*BetaVals(2)); % second Y derivative
            H2a=dip_convolve1d(toRegularize,[-1 0 1],0,1); %
            H12=dip_convolve1d(H2a,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
            myProjHessianSqr= H11 .*H11 + H22.*H22 + 2*H12.*H12; 
            T1 = 1+ alpha*(gamma*aReconSqr + myProjHessianSqr);  % will be summed over as log values
            tmp11 = dip_convolve1d(H11 ./ T1,s1*[1 -2 1],0,1)/(BetaVals(1)*BetaVals(1));
            tmp12=dip_convolve1d(H12 ./ T1,[-1 0 1],0,1); %
            tmp12=2*dip_convolve1d(tmp12,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
            tmp22 = dip_convolve1d((H22 ./ T1),s1*[1 -2 1],1,1)/(BetaVals(2)*BetaVals(2));
            myHessian2 = tmp11 + tmp12 + tmp22;
            myRegGrad = 2*alpha*myHessian2 + gamma*2*alpha*toRegularize ./ T1;
        else
            s1=1/1.6;s2=0.1727;
            H11=dip_convolve1d(toRegularize,s1*[1 -2 1],0,1)/(BetaVals(1)*BetaVals(1)); % second X derivative
            H22=dip_convolve1d(toRegularize,s1*[1 -2 1],1,1)/(BetaVals(2)*BetaVals(2)); % second Y derivative
            H33=dip_convolve1d(toRegularize,s1*[1 -2 1],2,1)/(BetaVals(3)*BetaVals(3)); % second Y derivative
            H2a=dip_convolve1d(toRegularize,[-1 0 1],0,1); %
            H12=dip_convolve1d(H2a,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
            H2b=dip_convolve1d(toRegularize,[-1 0 1],0,1); %
            H13=dip_convolve1d(H2b,s2*[-1 0 1],2,1)/(BetaVals(1)*BetaVals(3)); % mixed derivative XZ
            H2c=dip_convolve1d(toRegularize,[-1 0 1],1,1); %
            H23=dip_convolve1d(H2c,s2*[-1 0 1],2,1)/(BetaVals(2)*BetaVals(3)); % mixed derivative YZ
            myProjHessianSqr= H11 .*H11 + H22.*H22 +H33.*H33 + 2*H12.*H12+ 2*H13.*H13+ 2*H23.*H23;  % Does it need the weights of the filters?
            T1 = 1+ alpha*(gamma*aReconSqr + myProjHessianSqr);
            tmp11 = dip_convolve1d(H11 ./ T1,s1*[1 -2 1],0,1)/(BetaVals(1)*BetaVals(1));
            tmp22 = dip_convolve1d(H22 ./ T1,s1*[1 -2 1],1,1)/(BetaVals(2)*BetaVals(2));
            tmp33 = dip_convolve1d(H33 ./ T1,s1*[1 -2 1],1,1)/(BetaVals(3)*BetaVals(3));
            tmp12=dip_convolve1d(H12 ./ T1,[-1 0 1],0,1); %
            tmp12=2*dip_convolve1d(tmp12,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
            tmp13=dip_convolve1d(H13 ./ T1,[-1 0 1],0,1); %
            tmp13=2*dip_convolve1d(tmp13,s2*[-1 0 1],2,1)/(BetaVals(1)*BetaVals(3)); % mixed derivative XZ
            tmp23=dip_convolve1d(H23 ./ T1,[-1 0 1],1,1); %
            tmp23=2*dip_convolve1d(tmp23,s2*[-1 0 1],2,1)/(BetaVals(2)*BetaVals(3)); % mixed derivative YZ

            myHessian2 = tmp11 + tmp12 + tmp13 + tmp22+ tmp23+tmp33 ;
            myRegGrad = 2*alpha*myHessian2 + gamma*2*alpha*toRegularize ./ T1;
        end
        clear aReconSqr;
        % An alternative would have been to use the theorem:
        % (d/df)(d/dx)f(x) = f''(x) / f'(x)
        
        myReg = sum(log(T1));
