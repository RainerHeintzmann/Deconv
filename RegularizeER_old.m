% [myReg,myRegGrad]=RegularizeER(toRegularize,BetaVals) computes Arigovindan's roughness regularisation
% penalty=ln(epsR+ (|f|^2+ SumHessianSqr)),  SumHessianSqr= (d^2/dxdx f)^2 + (d^2/dydy f)^2+ (2d^2/dxdy f)^2
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient
% This function is based on a talk by
% Muthuvel  Arigovindan*,  Daniel Elnatan,  Jennifer  Fung,  Eric Branlund,  John W. Sedat,   and  David A.  Agard
% University of California at San Francisco
% suggesting a regularisation term as given below
%
% Now this is an assymetric implementation. The derivative is checked at the slightly wrong position

function [myReg,myRegGrad]=RegularizeER(toRegularize,BetaVals,epsR)
if nargin < 3
    epsR=1.0;
end
alpha= 1.0;  % Weighting between laplacian and noise penalty
%alpha= 10000.0;  % Weighting between laplacian and noise penalty
gamma = 1.0; % Weight for the intensity penalty under the logarithm
aReconSqr=abssqr(toRegularize);

if ndims(toRegularize) == 2 || (size(toRegularize,3) == 1)  % Two-D case
    % s1=1/1.6;s2=0.1727;   % These factors have been experimentally obtained by evaluating the rotation invariance of a test image
    s1=1;s2=sqrt(2);
    tRL1=circshift(toRegularize,[1 0]);tRR1=circshift(toRegularize,[-1 0]);
    tRL2=circshift(toRegularize,[0 1]);tRR2=circshift(toRegularize,[0 -1]);
    H11=s1*(tRL1+tRR1-2*toRegularize)/(BetaVals(1)*BetaVals(1)); % second X derivative
    H22=s1*(tRL2+tRR2-2*toRegularize)/(BetaVals(2)*BetaVals(2)); % second Y derivative
    if (0) % the gradient calculated below does not fit this case
        H12=s2*(-tRL1-tRL2+toRegularize+circshift(toRegularize,[-1 -1]))/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
    else
        tmp=(tRR1-tRL1)/(2*BetaVals(1));  % First order derivative d/dx f
        H12=(circshift(tmp,[0 -1])-circshift(tmp,[0 1]))/(2*BetaVals(2)); % mixed derivative XY
    end
    clear tRL1; clear tRL2; clear tRL3;
    myProjHessianSqr= abssqr(H11) + abssqr(H22) + 2*abssqr(H12);  % Calculate the sum of the squares of the invidual 2nd derivates
    T1 = epsR + alpha*(gamma*aReconSqr + myProjHessianSqr);
    % old gradient
%     tmp=H11;
%     tmp11 = s1*(circshift(tmp,[-1 0])-2*tmp+circshift(tmp,[1 0]))/(BetaVals(1)*BetaVals(1));
%     tmp=H22;
%     tmp22 = s1*(circshift(tmp,[0 -1])-2*tmp+circshift(tmp,[0 1]))/(BetaVals(2)*BetaVals(2));
%     if (0)
%         tmp=H12;
%         tmp12=s2*(-circshift(tmp,[-1 0])-circshift(tmp,[0 -1])+tmp+circshift(tmp,[-1 -1]))/(BetaVals(1)*BetaVals(2));
%     else
%         tmp=H12;
%         tmp=(circshift(tmp,[-1 0])-circshift(tmp,[1 0]))/(2*BetaVals(1));
%         tmp12=(circshift(tmp,[0 -1])-circshift(tmp,[0 1]))/(2*BetaVals(2));
%     end
%     myHessian2 = tmp11 + 2*tmp12 + tmp22;
%     %myHessian2 = tmp11 + tmp12 + tmp22;
%     myRegGrad = (2*alpha*myHessian2 + gamma*2*alpha*toRegularize) ./ T1;
    % Corrected gradient. s1 not in the gradient yet
    myRegGrad = alpha*((2*gamma*toRegularize - 4*H11/BetaVals(1)^2 - 4*H22/BetaVals(2)^2)./T1 + 2*(circshift(H11,[1 0])./circshift(T1, [1 0]) + ...
                circshift(H11, [-1 0])./circshift(T1, [-1 0]))/BetaVals(1)^2 + 2*(circshift(H22, [0 1])./circshift(T1, [0 1]) + ...
                circshift(H22, [0 -1])./circshift(T1, [0 -1]))/BetaVals(2)^2 + (circshift(H12, [1 1])./circshift(T1, [1 1]) + ...
                circshift(H12, [-1 -1])./circshift(T1, [-1 -1]) - circshift(H12, [1 -1])./circshift(T1, [1 -1]) - ...
                circshift(H12, [-1 1])./circshift(T1, [-1 1]))/BetaVals(1)/BetaVals(2));
else   % 3D case
    % s1=1/1.6;s2=0.1727;
    s1=1;s2=sqrt(2);
    tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
    tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
    tRL3=circshift(toRegularize,[0 0 1]);tRR3=circshift(toRegularize,[0 0 -1]);
    H11=s1*(tRL1+tRR1-2*toRegularize)/(BetaVals(1)*BetaVals(1)); % second X derivative
    H22=s1*(tRL2+tRR2-2*toRegularize)/(BetaVals(2)*BetaVals(2)); % second Y derivative
    H33=s1*(tRL3+tRR3-2*toRegularize)/(BetaVals(3)*BetaVals(3)); % second Z derivative
    if (0)  % assymmetric version; the gradient calculated below does not fit this case
        H12=s2*(-tRL1-tRL2+toRegularize+circshift(toRegularize,[1 1 0]))/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
        H13=s2*(-tRL1-tRL3+toRegularize+circshift(toRegularize,[1 0 1]))/(BetaVals(1)*BetaVals(3)); % mixed derivative XZ
        H23=s2*(-tRL2-tRL3+toRegularize+circshift(toRegularize,[0 1 1]))/(BetaVals(2)*BetaVals(3)); % mixed derivative YZ
    else % symmetric version
        tmp=(tRR1-tRL1)/(2*BetaVals(1));   % First order derivative d/dx f
        H12=(circshift(tmp,[0 -1 0])-circshift(tmp,[0 1 0]))/(2*BetaVals(2));
        H13=(circshift(tmp,[0 0 -1])-circshift(tmp,[0 0 1]))/(2*BetaVals(3));
        tmp=(tRR2-tRL2)/(2*BetaVals(2)); % First order derivative d/dy f
        H23=(circshift(tmp,[0 0 -1])-circshift(tmp,[0 0 1]))/(2*BetaVals(3));
    end
    clear tRL1; clear tRL2; clear tRL3; clear tmp;
    myProjHessianSqr= abssqr(H11) + abssqr(H22) + abssqr(H33) + 2*(abssqr(H12)+ abssqr(H13)+ abssqr(H23));  % Calculate the sum of the squares of the invidual 2nd derivates
    % myProjHessianSqr= H11 .*H11 + H22.*H22 +H33.*H33 + (H12.*H12+ H13.*H13+ H23.*H23);  % Does it need the weights of the filters?
    T1 = epsR + alpha*(gamma*aReconSqr + myProjHessianSqr);
%     if (0) % old gradient
%         tmp=H11;
%         tmp11 = s1*(circshift(tmp,[-1 0 0])-2*tmp+circshift(tmp,[1 0 0]))/(BetaVals(1)*BetaVals(1)); % d^2/(dxdx) H11/T1
%         tmp=H22;
%         tmp22 = s1*(circshift(tmp,[0 -1 0])-2*tmp+circshift(tmp,[0 1 0]))/(BetaVals(2)*BetaVals(2));
%         tmp=H33;
%         tmp33 = s1*(circshift(tmp,[0 0 -1])-2*tmp+circshift(tmp,[0 0 1]))/(BetaVals(3)*BetaVals(3));
%         if (0) % assymmetric version according to Arigovindan 2013
%             tmp=H12;
%             tmp12=s2*(-circshift(tmp,[1 0 0])-circshift(tmp,[0 1 0])+tmp+circshift(tmp,[1 1 0]))/(BetaVals(1)*BetaVals(2));
%             tmp=H13;
%             tmp13=s2*(-circshift(tmp,[1 0 0])-circshift(tmp,[0 0 1])+tmp+circshift(tmp,[1 0 1]))/(BetaVals(1)*BetaVals(3));
%             tmp=H23;
%             tmp23=s2*(-circshift(tmp,[0 1 0])-circshift(tmp,[0 0 1])+tmp+circshift(tmp,[0 1 1]))/(BetaVals(2)*BetaVals(3));
%         else  % symmetric version of second derivatives
%             tmp=H12;
%             tmp=(circshift(tmp,[-1 0 0])-circshift(tmp,[1 0 0]))/(2*BetaVals(1));
%             tmp12=(circshift(tmp,[0 -1 0])-circshift(tmp,[0 1 0]))/(2*BetaVals(2));
%             tmp=H13;
%             tmp=(circshift(tmp,[-1 0 0])-circshift(tmp,[1 0 0]))/(2*BetaVals(1));
%             tmp13=(circshift(tmp,[0 0 -1])-circshift(tmp,[0 0 1]))/(2*BetaVals(3));
%             tmp=H23;
%             tmp=(circshift(tmp,[0 -1 0])-circshift(tmp,[0 1 0]))/(2*BetaVals(2));
%             tmp23=(circshift(tmp,[0 0 -1])-circshift(tmp,[0 0 1]))/(2*BetaVals(3));
%         end
%         %             tmp12=dip_convolve1d(H12 ./ T1,[-1 0 1],0,1); %
%         %             tmp12=2*dip_convolve1d(tmp12,s2*[-1 0 1],1,1)/(BetaVals(1)*BetaVals(2)); % mixed derivative XY
%         %             tmp13=dip_convolve1d(H13 ./ T1,[-1 0 1],0,1); %
%         %             tmp13=2*dip_convolve1d(tmp13,s2*[-1 0 1],2,1)/(BetaVals(1)*BetaVals(3)); % mixed derivative XZ
%         %             tmp23=dip_convolve1d(H23 ./ T1,[-1 0 1],1,1); %
%         %             tmp23=2*dip_convolve1d(tmp23,s2*[-1 0 1],2,1)/(BetaVals(2)*BetaVals(3)); % mixed derivative YZ
%     
%         myHessian2 = tmp11 + tmp22 + tmp33 + 2*(tmp12 + tmp13 + tmp23);
%         % myHessian2 = tmp11 + tmp22 + tmp33 + (tmp12 + tmp13 + tmp23);
%         myRegGrad = (2*alpha*myHessian2 + gamma*2*alpha*toRegularize) ./ T1;
%         % myRegGrad = 2*alpha*myHessian2 + gamma*2*alpha*toRegularize ./ T1;
%     else

        % New gradient, not tested. Polina, 20.03.14
        myRegGrad = alpha*((2*gamma*toRegularize - 4*H11/BetaVals(1)^2 - 4*H22/BetaVals(2)^2 - 4*H33/BetaVals(3)^2)./T1 + ...
                    2*(circshift(H11,[1 0 0])./circshift(T1, [1 0 0]) + circshift(H11, [-1 0 0])./circshift(T1, [-1 0 0]))/BetaVals(1)^2 + ...
                    2*(circshift(H22, [0 1 0])./circshift(T1, [0 1 0]) + circshift(H22, [0 -1 0])./circshift(T1, [0 -1 0]))/BetaVals(2)^2 + ...
                    2*(circshift(H22, [0 0 1])./circshift(T1, [0 0 1]) + circshift(H22, [0 0 -1])./circshift(T1, [0 0 -1]))/BetaVals(3)^2 + ...
                    (circshift(H12, [1 1 0])./circshift(T1, [1 1 0]) + circshift(H12, [-1 -1 0])./circshift(T1, [-1 -1 0]) - ...
                    circshift(H12, [1 -1 0])./circshift(T1, [1 -1 0]) - circshift(H12, [-1 1 0])./circshift(T1, [-1 1 0]))/BetaVals(1)/BetaVals(2) + ...
                    (circshift(H13, [1 0 1])./circshift(T1, [1 0 1]) + circshift(H13, [-1 0 -1])./circshift(T1, [-1 0 -1]) - ...
                    circshift(H13, [1 0 -1])./circshift(T1, [1 0 -1]) - circshift(H13, [-1 0 1])./circshift(T1, [-1 0 1]))/BetaVals(1)/BetaVals(3) + ...
                    (circshift(H23, [0 1 1])./circshift(T1, [0 1 1]) + circshift(H23, [0 -1 -1])./circshift(T1, [0 -1 -1]) - ...
                    circshift(H23, [0 1 -1])./circshift(T1, [0 1 -1]) - circshift(H23, [0 -1 1])./circshift(T1, [0 -1 1]))/BetaVals(2)/BetaVals(3));
%     end
end
clear aReconSqr;
% An alternative would have been to use the theorem:
% (d/df)(d/dx)f(x) = f''(x) / f'(x)

myReg = sum(log(T1)); 
