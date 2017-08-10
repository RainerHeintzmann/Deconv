% [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR,doConvolve) computes Good's roughness regularisation:
% penalty= |Grad(f)|^2 / (|f|+epsR)
% toRegularize : 2D or 3D array to regularize
% BeatVals : vector of scaling factors (pixelsize)
% epsR : regularizes the division
% doConvolve : Flag that uses a convolved version of f in the nominator, if active
% myReg : Penalty value
% myRegGrad : Gradient
% See also: Verveer et al. Journal of Microscopy, 193, 50-61

% Corrected gradients, tested for 2D-case. Polina, 19.03.14

function [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR,doConvolve)
        % aGrad=gradient(aRecon)
        % myHessian = hessian(aRecon);
        % myRegGrad = -2*real(myHessian{1,1}+myHessian{2,2}+2*myHessian{1,2}); 
        if nargin < 3
            epsR=1e-4;
        end
        if nargin < 4
            doConvolve=0;   % If Active: This Trick devides by a better estimate of the local intensity.
        end
        if (ndims(toRegularize) == 2) || (size(toRegularize,3) == 1)
            if (ndims(toRegularize) == 2)
                tRL1=circshift(toRegularize,[1 0]);tRR1=circshift(toRegularize,[-1 0]);
                tRL2=circshift(toRegularize,[0 1]);tRR2=circshift(toRegularize,[0 -1]);
                aGrad{1}=(tRL1-tRR1)/(2*BetaVals(1));  % cyclic rotation
                aGrad{2}=(tRL2-tRR2)/(2*BetaVals(2));  % cyclic rotation
                if doConvolve
                    toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2)/5;  % blurs the nominator
                else
                    toRegularizeC=toRegularize;  
                end
                myRegGrad = (1/(2*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,[-2 0]))/(epsR + abs(tRR1)) + ...
                    (toRegularize - circshift(toRegularize,[2 0]))/(epsR + abs(tRL1))) + ...
                    (1/(2*BetaVals(2)*BetaVals(2)))*((toRegularize-circshift(toRegularize,[0 -2]))/(epsR + abs(tRR2)) + ...
                    (toRegularize - circshift(toRegularize,[0 2]))/(epsR + abs(tRL2))) - ...
                    sign(toRegularizeC)*(abssqr(aGrad{1})+abssqr(aGrad{2}))./((epsR+abs(toRegularizeC)).^2);
            else (size(toRegularize,3) == 1)
                tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
                tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
                if doConvolve
                    toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2)/5;  % blurs the nominator
                else
                    toRegularizeC=toRegularize;  
                end
                aGrad{1}=(tRL1-tRR1)/(2*BetaVals(1));  % cyclic rotation
                aGrad{2}=(tRL2-tRR2)/(2*BetaVals(2));  % cyclic rotation
                myRegGrad = (1/(2*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,[-2 0 0]))/(epsR + abs(tRR1)) + ...
                    (toRegularize - circshift(toRegularize,[2 0 0]))/(epsR + abs(tRL1))) + ...
                    (1/(2*BetaVals(2)*BetaVals(2)))*((toRegularize-circshift(toRegularize,[0 -2 0]))/(epsR + abs(tRR2)) + ...
                    (toRegularize - circshift(toRegularize,[0 2 0]))/(epsR + abs(tRL2))) - ...
                    sign(toRegularizeC)*(abssqr(aGrad{1})+abssqr(aGrad{2}))./((epsR+abs(toRegularizeC)).^2);
            end
        
            % myReg = sum(BetaVals(1)*aGrad{1} .* aGrad{1} + BetaVals(2)*aGrad{2} .* aGrad{2});
            myReg = sum((abssqr(aGrad{1}) + abssqr(aGrad{2})) ./ (epsR+abs(toRegularizeC)));
        elseif ndims(toRegularize) == 3
            tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
            tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
            tRL3=circshift(toRegularize,[0 0 1]);tRR3=circshift(toRegularize,[0 0 -1]);
            aGrad{1}=(tRL1-tRR1)/(2*BetaVals(1));  % cyclic rotation
            aGrad{2}=(tRL2-tRR2)/(2*BetaVals(2));  % cyclic rotation
            aGrad{3}=(tRL3-tRR3)/(2*BetaVals(3));  % cyclic rotation
            if doConvolve
                 toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
            else
                 toRegularizeC=toRegularize; 
            end
            myReg = sum((abssqr(aGrad{1}) + abssqr(aGrad{2}) + abssqr(aGrad{3})) ./ (epsR+abs(toRegularizeC)));
            myRegGrad = (1/(2*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,[-2 0 0])/(epsR + abs(tRR1))) + ...
                        (toRegularize-circshift(toRegularize,[2 0 0]))/(epsR + abs(tRL1))) + ...
                        (1/(2*BetaVals(2)*BetaVals(2)))*((toRegularize-circshift(toRegularize,[0 -2 0]))/(epsR + abs(tRR2)) + ...
                        (toRegularize-circshift(toRegularize,[0 2 0]))/(epsR + abs(tRL2))) + ...
                        (1/(2*BetaVals(3)*BetaVals(3)))*((toRegularize-circshift(toRegularize,[0 0 -2]))/(epsR + abs(tRR3)) + ...
                        (toRegularize-circshift(toRegularize,[0 0 2]))/(epsR + abs(tRL3))) - ...
                        sign(toRegularize)*(abssqr(aGrad{1})+abssqr(aGrad{2})+abssqr(aGrad{3}))./((epsR+abs(toRegularizeC)).^2);
        else % 1-D
            tRL1=circshift(toRegularize,1);tRR1=circshift(toRegularize,-1);
            aGrad{1}=(tRL1-tRR1)/(2*BetaVals(1));  % cyclic rotation
            if doConvolve
                 toRegularizeC=(toRegularize+tRL1+tRR1)/3;  % blurs the nominator
            else
                 toRegularizeC=toRegularize;  
            end
            myReg = sum(abssqr(aGrad{1}) ./ (epsR+abs(toRegularizeC)));
            myRegGrad = (1/(2*BetaVals(1)*BetaVals(1)))*((toRegularize-circshift(toRegularize,-2))/(epsR + abs(tRR1)) + ...
                        (toRegularize - circshift(toRegularize,2))/(epsR + abs(tRL1))) - ...
                        sign(toRegularize).*abssqr(aGrad{1})./((epsR+abs(toRegularizeC)).^2);
        end 
