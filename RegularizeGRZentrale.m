% [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR,doConvolve) computes Good's roughness regularisation:
% penalty= |Grad(f)|^2 / (|f|+epsR)
% toRegularize : 2D or 3D array to regularize
% BeatVals : vector of scaling factors (pixelsize)
% epsR : regularizes the division
% doConvolve : Flag that uses a convolved version of f in the nominator, if active.  DEPRECATED! Do not use!
% myReg : Penalty value
% myRegGrad : Gradient
% See also: Verveer et al. Journal of Microscopy, 193, 50-61

function [myReg,myRegGrad]=RegularizeGRZentrale(toRegularize,BetaVals,epsR,doConvolve)
        if nargin < 3
            epsR=0.1;  % If this value is too small the updates can lead to a numerical instability problem
        end
        if nargin < 4
            doConvolve=0;   % If Active: This Trick devides by a better estimate of the local intensity.
        end
        if doConvolve~=0
            error('doConvolve option is deprecated. Do not use.');
        end
        nom = epsR + abs(toRegularize);
        
tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
            tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
            tRL3=circshift(toRegularize,[0 0 1]);tRR3=circshift(toRegularize,[0 0 -1]);
            aGradZ{1}=(tRR1 - tRL1)/(2*BetaVals(1));
	    aGradZ{2}=(tRR2 - tRL2)/(2*BetaVals(2));
            aGradZ{3}=(tRR3 - tRL3)/(2*BetaVals(3));
%             if doConvolve
%                  toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
%             else
            if doConvolve
                 toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
            else
                 toRegularizeC=toRegularize; 
            end
           
            myReg = sum((abssqr(aGradZ{1}) + abssqr(aGradZ{2}) + abssqr(aGradZ{3})) ./...
                (epsR+abs(toRegularizeC)));
            
            nom = epsR + abs(toRegularizeC);
            myRegGrad = - ((abssqr(aGradZ{1}) + abssqr(aGradZ{2}) + abssqr(aGradZ{3})) ./...
                (nom.^2))+2*(circshift(aGradZ{1}, [1 0 0]))./circshift(nom, [1 0 0])-2*(circshift(aGradZ{1}, [-1 0 0]))./circshift(nom, [-1 0 0])+2*(circshift(aGradZ{2}, [0 1 0]))./circshift(nom, [0 1 0])-2*(circshift(aGradZ{2}, [0 -1 0]))./circshift(nom, [0 -1 0])+2*(circshift(aGradZ{3}, [0 0 1]))./circshift(nom, [0 0 1])-2*(circshift(aGradZ{3}, [0 0 -1]))./circshift(nom, [0 0 -1]);
       
  end
  %myRegGrad = 2*((aGradL{1}./circshift(nom, [1 0 0]) - aGradR{1}./circshift(nom, [-1 0 0]))/BetaVals(1) + ...
%                         (aGradL{2}./circshift(nom, [0 1 0]) - aGradR{2}./circshift(nom, [0 -1 0]))/BetaVals(2) + ...
%                         (aGradL{3}./circshift(nom, [0 0 1]) - aGradR{3}./circshift(nom, [0 0 -1]))/BetaVals(3) + ...
%                         ((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2) + (aGradL{3} - aGradR{3})/BetaVals(3))./nom) - ...
%                         sign(toRegularizeC).*(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2}) + ...
%                         abssqr(aGradL{3}) + abssqr(aGradR{3}))./nom.^2;