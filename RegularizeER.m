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
    if (ndims(toRegularize) == 2)
        tRL1=circshift(toRegularize,[1 0]);tRR1=circshift(toRegularize,[-1 0]);
        tRL2=circshift(toRegularize,[0 1]);tRR2=circshift(toRegularize,[0 -1]);
        tRL = circshift(toRegularize, [1 1]); tRR = circshift(toRegularize, [-1 -1]);
        tRLR = circshift(toRegularize, [1 -1]); tRRL = circshift(toRegularize, [-1 1]);
        H11=s1*(tRL1+tRR1-2*toRegularize)/(BetaVals(1)*BetaVals(1)); % second X derivative
        H22=s1*(tRL2+tRR2-2*toRegularize)/(BetaVals(2)*BetaVals(2)); % second Y derivative
        
        H12LL = (toRegularize + tRL - tRL1 - tRL2)/(BetaVals(1)*BetaVals(2));
        H12RR = (toRegularize + tRR - tRR1 - tRR2)/(BetaVals(1)*BetaVals(2));
        H12LR = (tRR2 + tRL1 - tRLR - toRegularize)/(BetaVals(1)*BetaVals(2));
        H12RL = (tRR1 + tRL2 - tRRL - toRegularize)/(BetaVals(1)*BetaVals(2));
    
        clear tRL1; clear tRL2;
        % Calculate the sum of the squares of the invidual 2nd derivates
        myProjHessianSqr= abssqr(H11) + abssqr(H22) + s2^2/4*(abssqr(H12LL) + abssqr(H12RR) + abssqr(H12LR) + abssqr(H12RL));  
        T1 = epsR + alpha*(gamma*aReconSqr + myProjHessianSqr);
        
        % gradient is tested
        myRegGrad = alpha*((2*gamma*toRegularize - 4*s1*H11/BetaVals(1)^2 - 4*s1*H22/BetaVals(2)^2)./T1 + ...
                    2*s1*(circshift(H11,[1 0])./circshift(T1, [1 0]) + circshift(H11, [-1 0])./circshift(T1, [-1 0]))/BetaVals(1)^2 + ...
                    2*s1*(circshift(H22, [0 1])./circshift(T1, [0 1]) + circshift(H22, [0 -1])./circshift(T1, [0 -1]))/BetaVals(2)^2 + ...
                    s2^2/2*((H12LL + H12RR - H12LR - H12RL)./T1 + circshift(H12RL - H12RR, [1 0])./circshift(T1, [1 0]) + ...
                    circshift(H12LR - H12LL, [-1 0])./circshift(T1, [-1 0]) + circshift(H12LR - H12RR, [0 1])./circshift(T1, [0 1]) + ...
                    circshift(H12RL - H12LL, [0 -1])./circshift(T1, [0 -1]) + circshift(H12LL, [-1 -1])./circshift(T1, [-1 -1]) + ...
                    circshift(H12RR, [1 1])./circshift(T1, [1 1]) - circshift(H12LR, [-1 1])./circshift(T1, [-1 1]) - ...
                    circshift(H12RL, [1 -1])./circshift(T1, [1 -1]))/BetaVals(1)/BetaVals(2));
    else
        tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
        tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
        tRL = circshift(toRegularize, [1 1 0]); tRR = circshift(toRegularize, [-1 -1 0]);
        tRLR = circshift(toRegularize, [1 -1 0]); tRRL = circshift(toRegularize, [-1 1 0]);
        H11=s1*(tRL1+tRR1-2*toRegularize)/(BetaVals(1)*BetaVals(1)); % second X derivative
        H22=s1*(tRL2+tRR2-2*toRegularize)/(BetaVals(2)*BetaVals(2)); % second Y derivative
        
        H12LL = (toRegularize + tRL - tRL1 - tRL2)/(BetaVals(1)*BetaVals(2));
        H12RR = (toRegularize + tRR - tRR1 - tRR2)/(BetaVals(1)*BetaVals(2));
        H12LR = (tRR2 + tRL1 - tRLR - toRegularize)/(BetaVals(1)*BetaVals(2));
        H12RL = (tRR1 + tRL2 - tRRL - toRegularize)/(BetaVals(1)*BetaVals(2));
    
        clear tRL1; clear tRL2;
        % Calculate the sum of the squares of the invidual 2nd derivates
        myProjHessianSqr= abssqr(H11) + abssqr(H22) + s2^2/4*(abssqr(H12LL) + abssqr(H12RR) + abssqr(H12LR) + abssqr(H12RL));  
        T1 = epsR + alpha*(gamma*aReconSqr + myProjHessianSqr);

        myRegGrad = alpha*((2*gamma*toRegularize - 4*s1*H11/BetaVals(1)^2 - 4*s1*H22/BetaVals(2)^2)./T1 + ...
                    2*s1*(circshift(H11,[1 0 0])./circshift(T1, [1 0 0]) + circshift(H11, [-1 0 0])./circshift(T1, [-1 0 0]))/BetaVals(1)^2 + ...
                    2*s1*(circshift(H22, [0 1 0])./circshift(T1, [0 1 0]) + circshift(H22, [0 -1 0])./circshift(T1, [0 -1 0]))/BetaVals(2)^2 + ...
                    s2^2/2*((H12LL + H12RR - H12LR - H12RL)./T1 + circshift(H12RL - H12RR, [1 0 0])./circshift(T1, [1 0 0]) + ...
                    circshift(H12LR - H12LL, [-1 0 0])./circshift(T1, [-1 0 0]) + circshift(H12LR - H12RR, [0 1 0])./circshift(T1, [0 1 0]) + ...
                    circshift(H12RL - H12LL, [0 -1 0])./circshift(T1, [0 -1 0]) + circshift(H12LL, [-1 -1 0])./circshift(T1, [-1 -1 0]) + ...
                    circshift(H12RR, [1 1 0])./circshift(T1, [1 1 0]) - circshift(H12LR, [-1 1 0])./circshift(T1, [-1 1 0]) - ...
                    circshift(H12RL, [1 -1 0])./circshift(T1, [1 -1 0]))/BetaVals(1)/BetaVals(2));
    end
else   % 3D case
    % s1=1/1.6;s2=0.1727;
    s1=1;s2=sqrt(2);
    tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
    tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
    tRL3=circshift(toRegularize,[0 0 1]);tRR3=circshift(toRegularize,[0 0 -1]);
    tRL12 = circshift(toRegularize, [1 1 0]); tRR12 = circshift(toRegularize, [-1 -1 0]);
    tRLR12 = circshift(toRegularize, [1 -1 0]); tRRL12 = circshift(toRegularize, [-1 1 0]);
    tRL13 = circshift(toRegularize, [1 0 1]); tRR13 = circshift(toRegularize, [-1 0 -1]);
    tRLR13 = circshift(toRegularize, [1 0 -1]); tRRL13 = circshift(toRegularize, [-1 0 1]);
    tRL23 = circshift(toRegularize, [0 1 1]); tRR23 = circshift(toRegularize, [0 -1 -1]);
    tRLR23 = circshift(toRegularize, [0 1 -1]); tRRL23 = circshift(toRegularize, [0 -1 1]);
    H11=s1*(tRL1+tRR1-2*toRegularize)/(BetaVals(1)*BetaVals(1)); % second X derivative
    H22=s1*(tRL2+tRR2-2*toRegularize)/(BetaVals(2)*BetaVals(2)); % second Y derivative
    H33=s1*(tRL3+tRR3-2*toRegularize)/(BetaVals(3)*BetaVals(3)); % second Y derivative
    
    H12LL = (toRegularize + tRL12 - tRL1 - tRL2)/(BetaVals(1)*BetaVals(2));
    H12RR = (toRegularize + tRR12 - tRR1 - tRR2)/(BetaVals(1)*BetaVals(2));
    H12LR = (tRR2 + tRL1 - tRLR12 - toRegularize)/(BetaVals(1)*BetaVals(2));
    H12RL = (tRR1 + tRL2 - tRRL12 - toRegularize)/(BetaVals(1)*BetaVals(2));
    
    H13LL = (toRegularize + tRL13 - tRL1 - tRL3)/(BetaVals(1)*BetaVals(3));
    H13RR = (toRegularize + tRR13 - tRR1 - tRR3)/(BetaVals(1)*BetaVals(3));
    H13LR = (tRR3 + tRL1 - tRLR13 - toRegularize)/(BetaVals(1)*BetaVals(3));
    H13RL = (tRR1 + tRL3 - tRRL13 - toRegularize)/(BetaVals(1)*BetaVals(3));

    H23LL = (toRegularize + tRL23 - tRL2 - tRL3)/(BetaVals(2)*BetaVals(3));
    H23RR = (toRegularize + tRR23 - tRR2 - tRR3)/(BetaVals(2)*BetaVals(3));
    H23LR = (tRR3 + tRL2 - tRLR23 - toRegularize)/(BetaVals(2)*BetaVals(3));
    H23RL = (tRR2 + tRL3 - tRRL23 - toRegularize)/(BetaVals(2)*BetaVals(3));
    
    clear tRL1; clear tRL2; clear tRL3; clear tmp;
    % Calculate the sum of the squares of the invidual 2nd derivates
    myProjHessianSqr =  abssqr(H11) + abssqr(H22) + abssqr(H33) + s2^2/4*(abssqr(H12LL) + abssqr(H12RR) + abssqr(H12LR) + abssqr(H12RL) + ...
                        abssqr(H13LL) + abssqr(H13RR) + abssqr(H13LR) + abssqr(H13RL) + abssqr(H23LL) + abssqr(H23RR) + abssqr(H23LR) + abssqr(H23RL));  
    % myProjHessianSqr= H11 .*H11 + H22.*H22 +H33.*H33 + (H12.*H12+ H13.*H13+ H23.*H23);  % Does it need the weights of the filters?
    T1 = epsR + alpha*(gamma*aReconSqr + myProjHessianSqr);
    
    myRegGrad = alpha*((2*gamma*toRegularize - (4*s1)*H11/BetaVals(1)^2 - (4*s1)*H22/BetaVals(2)^2 - (4*s1)*H33/BetaVals(3)^2)./T1 + ...
                (2*s1)*(circshift(H11,[1 0 0])./circshift(T1, [1 0 0]) + circshift(H11, [-1 0 0])./circshift(T1, [-1 0 0]))/BetaVals(1)^2 + ...
                (2*s1)*(circshift(H22, [0 1 0])./circshift(T1, [0 1 0]) + circshift(H22, [0 -1 0])./circshift(T1, [0 -1 0]))/BetaVals(2)^2 + ...
                (2*s1)*(circshift(H33, [0 0 1])./circshift(T1, [0 0 1]) + circshift(H33, [0 0 -1])./circshift(T1, [0 0 -1]))/BetaVals(3)^2 + ...
                (s2^2/2)*((H12LL + H12RR - H12LR - H12RL)./T1 + circshift(H12RL - H12RR, [1 0 0])./circshift(T1, [1 0 0]) + ...
                circshift(H12LR - H12LL, [-1 0 0])./circshift(T1, [-1 0 0]) + circshift(H12LR - H12RR, [0 1 0])./circshift(T1, [0 1 0]) + ...
                circshift(H12RL - H12LL, [0 -1 0])./circshift(T1, [0 -1 0]) + circshift(H12LL, [-1 -1 0])./circshift(T1, [-1 -1 0]) + ...
                circshift(H12RR, [1 1 0])./circshift(T1, [1 1 0]) - circshift(H12LR, [-1 1 0])./circshift(T1, [-1 1 0]) - ...
                circshift(H12RL, [1 -1 0])./circshift(T1, [1 -1 0]))/BetaVals(1)/BetaVals(2) + ...
                (s2^2/2)*((H13LL + H13RR - H13LR - H13RL)./T1 + circshift(H13RL - H13RR, [1 0 0])./circshift(T1, [1 0 0]) + ...
                circshift(H13LR - H13LL, [-1 0 0])./circshift(T1, [-1 0 0]) + circshift(H13LR - H13RR, [0 0 1])./circshift(T1, [0 0 1]) + ...
                circshift(H13RL - H13LL, [0 0 -1])./circshift(T1, [0 0 -1]) + circshift(H13LL, [-1 0 -1])./circshift(T1, [-1 0 -1]) + ...
                circshift(H13RR, [1 0 1])./circshift(T1, [1 0 1]) - circshift(H13LR, [-1 0 1])./circshift(T1, [-1 0 1]) - ...
                circshift(H13RL, [1 0 -1])./circshift(T1, [1 0 -1]))/BetaVals(1)/BetaVals(3) + ...
                (s2^2/2)*((H23LL + H23RR - H23LR - H23RL)./T1 + circshift(H23RL - H23RR, [0 1 0])./circshift(T1, [0 1 0]) + ...
                circshift(H23LR - H23LL, [0 -1 0])./circshift(T1, [0 -1 0]) + circshift(H23LR - H23RR, [0 0 1])./circshift(T1, [0 0 1]) + ...
                circshift(H23RL - H23LL, [0 0 -1])./circshift(T1, [0 0 -1]) + circshift(H23LL, [0 -1 -1])./circshift(T1, [0 -1 -1]) + ...
                circshift(H23RR, [0 1 1])./circshift(T1, [0 1 1]) - circshift(H23LR, [0 -1 1])./circshift(T1, [0 -1 1]) - ...
                circshift(H23RL, [0 1 -1])./circshift(T1, [0 1 -1]))/BetaVals(2)/BetaVals(3));
    end
clear aReconSqr;
% An alternative would have been to use the theorem:
% (d/df)(d/dx)f(x) = f''(x) / f'(x)

myReg = sum(log(T1)); 
