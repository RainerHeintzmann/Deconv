% [myReg,myRegGrad]=RegularizeGRLap6(toRegularize,BetaVals,epsR,doConvolve) computes Good's roughness regularisation:   A 6-point Laplacian regularizer with Good's roughness scaling
% penalty= |Grad(f)|^2 / (|f|+epsR)
% toRegularize : 2D or 3D array to regularize
% BeatVals : vector of scaling factors (pixelsize)
% epsR : regularizes the division
% doConvolve : Flag that uses a convolved version of f in the nominator, if active.  DEPRECATED! Do not use!
% myReg : Penalty value
% myRegGrad : Gradient
% See also: Verveer et al. Journal of Microscopy, 193, 50-61
%
%  Author of this regularizer: Felix Mucha
%

function [myReg,myRegGrad]=RegularizeGRLap6(toRegularize,BetaVals,epsR,doConvolve)
        if nargin < 3
            epsR=0.1;  % If this value is too small the updates can lead to a numerical instability problem
        end
        if nargin < 4
            doConvolve=0;   % If Active: This Trick devides by a better estimate of the local intensity.
        end
        if doConvolve~=0
            error('doConvolve option is deprecated. Do not use.');
        end
        epsR=0.1;
        nom = epsR + abs(toRegularize);
        
tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
            tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
            tRL3=circshift(toRegularize,[0 0 1]);tRR3=circshift(toRegularize,[0 0 -1]);
            aGradLap{1}=(tRL1+tRR1-2*toRegularize)/(BetaVals(1)).^2;
            aGradLap{2}=(tRL2+tRR2-2*toRegularize)/(BetaVals(2)).^2;
            aGradLap{3}=(tRL3+tRR3-2*toRegularize)/(BetaVals(3)).^2;
%             aGradZ{1}=(tRR1 - tRL1)/(2*BetaVals(1));
%             
% 	    aGradZ{2}=(tRR2 - tRL2)/(2*BetaVals(2));
%             aGradZ{3}=(tRR3 - tRL3)/(2*BetaVals(3));
%             if doConvolve
%                  toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
%             else
            if doConvolve
                 toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
            else
                 toRegularizeC=toRegularize; 
            end
            toRegularizeC=toRegularize;
            nom = (abs(toRegularizeC+epsR));
            myReg = sum(((aGradLap{1}+aGradLap{2}+aGradLap{3}))./nom);
            
            myRegGrad1X=0;
            myRegGrad2X=0;
            myRegGrad3X=0;
            myRegGrad1Y=0;
            myRegGrad2Y=0;
            myRegGrad3Y=0;
            myRegGrad1Z=0;
            myRegGrad2Z=0;
            myRegGrad3Z=0;
           
            myRegGrad1X=myRegGrad1X+(-2)./nom./(BetaVals(1)).^2;
            myRegGrad1Y=myRegGrad1Y+(-2)./nom./(BetaVals(2)).^2;
            myRegGrad1Z=myRegGrad1Z+(-2)./nom./(BetaVals(3)).^2;
            myRegGrad2X=myRegGrad2X+1./(circshift(nom,[-1 0 0]))./(BetaVals(1)).^2+1./(circshift(nom,[1 0 0]))./(BetaVals(1)).^2;
            myRegGrad2Y=myRegGrad2Y+1./(circshift(nom,[0 -1 0]))./(BetaVals(2)).^2+1./(circshift(nom,[0 1 0]))./(BetaVals(2)).^2;
            myRegGrad2Z=myRegGrad2Z+1./circshift(nom,[0 0 -1])./(BetaVals(3)).^2+1./(circshift(nom,[0 0 1]))./(BetaVals(3)).^2;
            myRegGrad3X= myRegGrad3X-sign(nom).*((aGradLap{1}))./nom;
            myRegGrad3Y= myRegGrad3Y-sign(nom).*((aGradLap{2}))./nom;
            myRegGrad3Z= myRegGrad3Z-sign(nom).*((aGradLap{3}))./nom;

            myRegGrad=myRegGrad2X+myRegGrad2Y+myRegGrad2Z+myRegGrad1X+myRegGrad1Y+myRegGrad1Z+(myRegGrad3X+myRegGrad3Y+myRegGrad3Z)./nom;  

   end           


  %myRegGrad = 2*((aGradL{1}./circshift(nom, [1 0 0]) - aGradR{1}./circshift(nom, [-1 0 0]))/BetaVals(1) + ...
%                         (aGradL{2}./circshift(nom, [0 1 0]) - aGradR{2}./circshift(nom, [0 -1 0]))/BetaVals(2) + ...
%                         (aGradL{3}./circshift(nom, [0 0 1]) - aGradR{3}./circshift(nom, [0 0 -1]))/BetaVals(3) + ...
%                         ((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2) + (aGradL{3} - aGradR{3})/BetaVals(3))./nom) - ...
%                         sign(toRegularizeC).*(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2}) + ...
%                         abssqr(aGradL{3}) + abssqr(aGradR{3}))./nom.^2;