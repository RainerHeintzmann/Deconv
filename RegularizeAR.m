% [myReg,myRegGrad]=RegularizeAR(toRegularize,BetaVals) computes Arigovindan's roughness regularisation
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient
% This function is based on a talk by 
% Muthuvel  Arigovindan*,  Daniel Elnatan,  Jennifer  Fung,  Eric Branlund,  John W. Sedat,   and  David A.  Agard
% University of California at San Francisco
% suggesting a regolarisation term as given below
%
% Now this is an assymetric implementation. The derivative is checked at the slightly wrong position

function [myReg,myRegGrad]=RegularizeAR(toRegularize,BetaVals,epsR)
if nargin < 3
    epsR=1.0;
end
        alpha= 1.0;  % Weighting between laplacian and noise penalty
        %alpha= 10000.0;  % Weighting between laplacian and noise penalty
        gamma = 1.0; % Weight for the intensity penalty under the logarithm
        aReconSqr=toRegularize.*toRegularize;

        if ndims(toRegularize) == 2 || (size(toRegularize,3) == 1)
            % s1=1/1.6;s2=0.1727;
            s1=1;s2=sqrt(2);
            tRL1=circshift(toRegularize,[1 0]);tRR1=circshift(toRegularize,[-1 0]); 
            tRL2=circshift(toRegularize,[0 1]); 
            H11=s1*(tRL1+tRR1-2*toRegularize)/(BetaVals(1)*BetaVals(1)); % second X derivative
            H22=s1*(tRL2+circshift(toRegularize,[0 -1])-2*toRegularize)/(BetaVals(2)*BetaVals(2)); % second Y derivative
            if (0)
                H12=s2*(-tRL1-tRL2+toRegularize+circshift(toRegularize,[-1 -1]))/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
            else
                tmp=(tRR1-tRL1)/(2*BetaVals(1));
                H12=(circshift(tmp,[0 -1])-circshift(tmp,[0 1]))/(2*BetaVals(2));
            end
            clear tRL1; clear tRL2; clear tRL3;
            myProjHessianSqr= H11 .*H11 + H22.*H22 + 2*H12.*H12;  % Does it need the weights of the filters?
            %myProjHessianSqr= H11 .*H11 + H22.*H22 + H12.*H12;  % Does it need the weights of the filters?
            T1 = epsR + alpha*(gamma*aReconSqr + myProjHessianSqr);
            tmp=H11./T1;
            tmp11 = s1*(circshift(tmp,[-1 0])-2*tmp+circshift(tmp,[1 0]))/(BetaVals(1)*BetaVals(1));
            tmp=H22./T1;
            tmp22 = s1*(circshift(tmp,[0 -1])-2*tmp+circshift(tmp,[0 1]))/(BetaVals(2)*BetaVals(2));
            if (0)
                tmp=H12./T1;
                tmp12=s2*(-circshift(tmp,[-1 0])-circshift(tmp,[0 -1])+tmp+circshift(tmp,[-1 -1]))/(BetaVals(1)*BetaVals(2));
            else
                tmp=H12./T1;
                tmp=(circshift(tmp,[-1 0])-circshift(tmp,[1 0]))/(2*BetaVals(1));
                tmp12=(circshift(tmp,[0 -1])-circshift(tmp,[0 1]))/(2*BetaVals(2));
            end
            myHessian2 = tmp11 + 2*tmp12 + tmp22;
            %myHessian2 = tmp11 + tmp12 + tmp22;
            myRegGrad = (2*alpha*myHessian2 + gamma*2*alpha*toRegularize) ./ T1;
        else
            % s1=1/1.6;s2=0.1727;
            s1=1;s2=sqrt(2);
            tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]); 
            tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]); 
            tRL3=circshift(toRegularize,[0 0 1]); 
            H11=s1*(tRL1+tRR1-2*toRegularize)/(BetaVals(1)*BetaVals(1)); % second X derivative
            H22=s1*(tRL2+tRR2-2*toRegularize)/(BetaVals(2)*BetaVals(2)); % second Y derivative
            H33=s1*(tRL3+circshift(toRegularize,[0 0 -1])-2*toRegularize)/(BetaVals(3)*BetaVals(3)); % second Y derivative           
            if (0)  % assymmetric version
                H12=s2*(-tRL1-tRL2+toRegularize+circshift(toRegularize,[1 1 0]))/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
                H13=s2*(-tRL1-tRL3+toRegularize+circshift(toRegularize,[1 0 1]))/(BetaVals(1)*BetaVals(3)); % mixed derivative XZ
                H23=s2*(-tRL2-tRL3+toRegularize+circshift(toRegularize,[0 1 1]))/(BetaVals(2)*BetaVals(3)); % mixed derivative YZ
            else % symmetric version
                tmp=(tRR1-tRL1)/(2*BetaVals(1));
                H12=(circshift(tmp,[0 -1 0])-circshift(tmp,[0 1 0]))/(2*BetaVals(2));
                H13=(circshift(tmp,[0 0 -1])-circshift(tmp,[0 0 1]))/(2*BetaVals(3));
                H23=(tRR2-tRL2)/(2*BetaVals(2));
                H23=(circshift(H23,[0 0 -1])-circshift(H23,[0 0 1]))/(2*BetaVals(3));
            end
            clear tRL1; clear tRL2; clear tRL3;
            myProjHessianSqr= H11 .*H11 + H22.*H22 +H33.*H33 + 2*(H12.*H12+ H13.*H13+ H23.*H23);  % Does it need the weights of the filters?
            % myProjHessianSqr= H11 .*H11 + H22.*H22 +H33.*H33 + (H12.*H12+ H13.*H13+ H23.*H23);  % Does it need the weights of the filters?
            T1 = epsR + alpha*(gamma*aReconSqr + myProjHessianSqr);
            tmp=H11./T1;
            tmp11 = s1*(circshift(tmp,[-1 0 0])-2*tmp+circshift(tmp,[1 0 0]))/(BetaVals(1)*BetaVals(1));
            tmp=H22./T1;
            tmp22 = s1*(circshift(tmp,[0 -1 0])-2*tmp+circshift(tmp,[0 1 0]))/(BetaVals(2)*BetaVals(2));
            tmp=H33./T1;
            tmp33 = s1*(circshift(tmp,[0 0 -1])-2*tmp+circshift(tmp,[0 0 1]))/(BetaVals(3)*BetaVals(3));
            if (0) % assymmetric version according to Arigovindan 2013
                tmp=H12./T1;
                tmp12=s2*(-circshift(tmp,[1 0 0])-circshift(tmp,[0 1 0])+tmp+circshift(tmp,[1 1 0]))/(BetaVals(1)*BetaVals(2));
                tmp=H13./T1;
                tmp13=s2*(-circshift(tmp,[1 0 0])-circshift(tmp,[0 0 1])+tmp+circshift(tmp,[1 0 1]))/(BetaVals(1)*BetaVals(3));
                tmp=H23./T1;
                tmp23=s2*(-circshift(tmp,[0 1 0])-circshift(tmp,[0 0 1])+tmp+circshift(tmp,[0 1 1]))/(BetaVals(2)*BetaVals(3));
            else  % symmetric version of second derivatives
                tmp=H12./T1;
                tmp=(circshift(tmp,[-1 0 0])-circshift(tmp,[1 0 0]))/(2*BetaVals(1));
                tmp12=(circshift(tmp,[0 -1 0])-circshift(tmp,[0 1 0]))/(2*BetaVals(2));
                tmp=H13./T1;
                tmp=(circshift(tmp,[-1 0 0])-circshift(tmp,[1 0 0]))/(2*BetaVals(1));
                tmp13=(circshift(tmp,[0 0 -1])-circshift(tmp,[0 0 1]))/(2*BetaVals(3));
                tmp=H23./T1;
                tmp=(circshift(tmp,[0 -1 0])-circshift(tmp,[0 1 0]))/(2*BetaVals(2));
                tmp23=(circshift(tmp,[0 0 -1])-circshift(tmp,[0 0 1]))/(2*BetaVals(3));
            end
%             tmp12=dip_convolve1d(H12 ./ T1,[-1 0 1],0,1); %
%             tmp12=2*dip_convolve1d(tmp12,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
%             tmp13=dip_convolve1d(H13 ./ T1,[-1 0 1],0,1); %
%             tmp13=2*dip_convolve1d(tmp13,s2*[-1 0 1],2,1)/(BetaVals(1)*BetaVals(3)); % mixed derivative XZ
%             tmp23=dip_convolve1d(H23 ./ T1,[-1 0 1],1,1); %
%             tmp23=2*dip_convolve1d(tmp23,s2*[-1 0 1],2,1)/(BetaVals(2)*BetaVals(3)); % mixed derivative YZ

            myHessian2 = tmp11 + tmp22 + tmp33 + 2*(tmp12 + tmp13 + tmp23);
             % myHessian2 = tmp11 + tmp22 + tmp33 + (tmp12 + tmp13 + tmp23);
            myRegGrad = (2*alpha*myHessian2 + gamma*2*alpha*toRegularize) ./ T1;
            % myRegGrad = 2*alpha*myHessian2 + gamma*2*alpha*toRegularize ./ T1;
        end
        clear aReconSqr;
        % An alternative would have been to use the theorem:
        % (d/df)(d/dx)f(x) = f''(x) / f'(x)
        
        myReg = sum(log(T1));
