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
            
    	    aGradX=(tRR1 - tRL1)./(2*BetaVals(1));
            aGradY=(tRR2 - tRL2)./(2*BetaVals(2));
            aGradZ=(tRR3 - tRL3)./(2*BetaVals(3));
%             aGradZ{3}=(tRR3 - tRL3)/(BetaVals(3)).^2;
%             if doConvolve
%                  toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
%             else
            if doConvolve
                 toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
            else
                 toRegularizeC=toRegularize; 
            end
            myRegGrad1X=0;
            myRegGrad2X=0;
            myRegGrad3X=0;
            myRegGrad1Y=0;
            myRegGrad2Y=0;
            myRegGrad3Y=0;
            myRegGrad1Z=0;
            myRegGrad2Z=0;
            myRegGrad3Z=0;

            nom = epsR + abs(toRegularizeC);
            myReg = sum((((abssqr(aGradX)+abssqr(aGradY)+abssqr(aGradZ))./nom))  );
                       
            myRegGrad1X =myRegGrad1X -sign(nom).* (abssqr(aGradX))./(nom);
            myRegGrad1Y =myRegGrad1Y -sign(nom).* (abssqr(aGradY))./(nom);
            myRegGrad1Z =myRegGrad1Z -sign(nom).* (abssqr(aGradZ))./(nom);
            myRegGrad2X=myRegGrad2X+(((circshift(aGradX, [1 0 0]))))./(2*BetaVals(1))./(circshift(nom, [1 0 0]));
            myRegGrad2Y=myRegGrad2Y+(((circshift(aGradY, [0 1 0]))))./(2*BetaVals(2))./(circshift(nom, [0 1 0]));
            myRegGrad2Z=myRegGrad2Z+(((circshift(aGradZ, [0 0 1]))))./(2*BetaVals(3))./(circshift(nom, [0 0 1]));
            myRegGrad3X= myRegGrad3X-(((circshift(aGradX, [-1 0 0]))))./(2*BetaVals(1))./(circshift(nom, [-1 0 0]));
            myRegGrad3Y= myRegGrad3Y-(((circshift(aGradY, [0 -1 0]))))./(2*BetaVals(2))./(circshift(nom, [0 -1 0]));
            myRegGrad3Z= myRegGrad3Z-(((circshift(aGradZ, [0 0 -1]))))./(2*BetaVals(3))./(circshift(nom, [0 0 -1]));
myRegGrad=2*myRegGrad2X+2*myRegGrad3X+2*myRegGrad2Y+2*myRegGrad3Y+2*myRegGrad2Z+2*myRegGrad3Z+(myRegGrad1X+myRegGrad1Y+myRegGrad1Z)./nom;  
end
  %myRegGrad = 2*((aGradL{1}./circshift(nom, [1 0 0]) - aGradR{1}./circshift(nom, [-1 0 0]))/BetaVals(1) + ...
%                         (aGradL{2}./circshift(nom, [0 1 0]) - aGradR{2}./circshift(nom, [0 -1 0]))/BetaVals(2) + ...
%                         (aGradL{3}./circshift(nom, [0 0 1]) - aGradR{3}./circshift(nom, [0 0 -1]))/BetaVals(3) + ...
%                         ((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2) + (aGradL{3} - aGradR{3})/BetaVals(3))./nom) - ...
%                         sign(toRegularizeC).*(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2}) + ...
%                         abssqr(aGradL{3}) + abssqr(aGradR{3}))./nom.^2;