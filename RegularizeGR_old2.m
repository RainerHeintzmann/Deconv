% [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR,doConvolve) computes Good's roughness regularisation:
% penalty= |Grad(f)|^2 / (|f|+epsR)
% toRegularize : 2D or 3D array to regularize
% BeatVals : vector of scaling factors (pixelsize)
% epsR : regularizes the division
% doConvolve : Flag that uses a convolved version of f in the nominator, if active
% myReg : Penalty value
% myRegGrad : Gradient
% See also: Verveer et al. Journal of Microscopy, 193, 50-61

function [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR,doConvolve)
        % aGrad=gradient(aRecon)
        % myHessian = hessian(aRecon);
        % myRegGrad = -2*real(myHessian{1,1}+myHessian{2,2}+2*myHessian{1,2}); 
        if nargin < 3
            epsR=0.1;  % If this value is too small the updates can lead to a numerical instability problem
        end
        if nargin < 4
            doConvolve=0;   % If Active: This Trick devides by a better estimate of the local intensity.
        end
        if (ndims(toRegularize) == 2) || (size(toRegularize,3) == 1)
            if (ndims(toRegularize) == 2)
                tRL1=circshift(toRegularize,[1 0]);tRR1=circshift(toRegularize,[-1 0]);
                tRL2=circshift(toRegularize,[0 1]);tRR2=circshift(toRegularize,[0 -1]);
                aGradL{1}=(toRegularize - tRL1)/BetaVals(1);  % cyclic rotation
                aGradL{2}=(toRegularize - tRL2)/BetaVals(2);  % cyclic rotation
                aGradR{1}=(tRR1 - toRegularize)/BetaVals(1);  % cyclic rotation
                aGradR{2}=(tRR2 - toRegularize)/BetaVals(2);  % cyclic rotation
                if doConvolve
                    toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2)/5;  % blurs the nominator
                else
                    toRegularizeC=toRegularize;  
                end
                nom = epsR + abs(toRegularizeC);
                myRegGrad = 2*((aGradL{1}./circshift(nom, [1 0]) - aGradR{1}./circshift(nom, [-1 0]))/BetaVals(1) + ...
                            (aGradL{2}./circshift(nom, [0 1]) - aGradR{2}./circshift(nom, [0 -1]))/BetaVals(2) + ...
                            ((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2))./nom) - ...
                            sign(toRegularizeC).*(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2}))./nom.^2;
            else (size(toRegularize,3) == 1)
                tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
                tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
                if doConvolve
                    toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2)/5;  % blurs the nominator
                else
                    toRegularizeC=toRegularize;  
                end
                aGradL{1}=(toRegularize - tRL1)/BetaVals(1);  % cyclic rotation
                aGradL{2}=(toRegularize - tRL2)/BetaVals(2);  % cyclic rotation
                aGradR{1}=(tRR1 - toRegularize)/BetaVals(1);  % cyclic rotation
                aGradR{2}=(tRR2 - toRegularize)/BetaVals(2);  % cyclic rotation 
                nom = epsR + abs(toRegularizeC);
                myRegGrad = 2*((aGradL{1}./circshift(nom, [1 0 0]) - aGradR{1}./circshift(nom, [-1 0 0]))/BetaVals(1) + ...
                            (aGradL{2}./circshift(nom, [0 1 0]) - aGradR{2}./circshift(nom, [0 -1 0]))/BetaVals(2) + ...
                            ((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2))./nom) - ...
                            sign(toRegularizeC).*(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2}))./nom.^2;
            end
            % myReg = sum(BetaVals(1)*aGrad{1} .* aGrad{1} + BetaVals(2)*aGrad{2} .* aGrad{2});
            myReg = sum((abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2})) ./ (epsR+abs(toRegularizeC)));
        elseif ndims(toRegularize) == 3
            tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
            tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
            tRL3=circshift(toRegularize,[0 0 1]);tRR3=circshift(toRegularize,[0 0 -1]);
            aGradL{1}=(toRegularize - tRL1)/BetaVals(1);  % cyclic rotation
            aGradL{2}=(toRegularize - tRL2)/BetaVals(2);  % cyclic rotation
            aGradL{3}=(toRegularize - tRL3)/BetaVals(3);  % cyclic rotation
            aGradR{1}=(tRR1 - toRegularize)/BetaVals(1);  % cyclic rotation
            aGradR{2}=(tRR2 - toRegularize)/BetaVals(2);  % cyclic rotation
            aGradR{3}=(tRR3 - toRegularize)/BetaVals(3);  % cyclic rotation
            if doConvolve
                 toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
            else
                 toRegularizeC=toRegularize; 
            end
            myReg = sum((abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradL{3}) + abssqr(aGradR{1}) + abssqr(aGradR{2}) + abssqr(aGradR{3})) ./...
                (epsR+abs(toRegularizeC)));
            nom = epsR + abs(toRegularizeC);
            myRegGrad = 2*((aGradL{1}./circshift(nom, [1 0 0]) - aGradR{1}./circshift(nom, [-1 0 0]))/BetaVals(1) + ...
                        (aGradL{2}./circshift(nom, [0 1 0]) - aGradR{2}./circshift(nom, [0 -1 0]))/BetaVals(2) + ...
                        (aGradL{3}./circshift(nom, [0 0 1]) - aGradR{3}./circshift(nom, [0 0 -1]))/BetaVals(3) + ...
                        ((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2) + (aGradL{3} - aGradR{3})/BetaVals(3))./nom) - ...
                        sign(toRegularizeC).*(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2}) + ...
                        abssqr(aGradL{3}) + abssqr(aGradR{3}))./nom.^2;
        else % 1-D
            tRL1=circshift(toRegularize,1);tRR1=circshift(toRegularize,-1);
            aGradL=(toRegularize - tRL1)/BetaVals(1);  % cyclic rotation
            aGradR=(tRR1 - toRegularize)/BetaVals(1);  % cyclic rotation
            if doConvolve
                 toRegularizeC=(toRegularize+tRL1+tRR1)/3;  % blurs the nominator
            else
                 toRegularizeC=toRegularize;  
            end
            myReg = sum((abssqr(aGradL) + abssqr(aGradR)) ./ (epsR+abs(toRegularizeC)));
            nom = epsR + abs(toRegularizeC);
            myRegGrad = 2*(aGradL./circshift(nom, 1) - aGradR./circshift(nom, -1)  + (aGradL - aGradR)./nom)/BetaVals(1)...
                        - sign(toRegularizeC).*(abssqr(aGradL) + abssqr(aGradR))./nom.^2;
        end 
