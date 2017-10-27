% [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR,doConvolve) computes Good's roughness regularisation:
% penalty= |Grad(f)|^2 / (|f|+epsR)
% toRegularize : 2D or 3D array to regularize
% BeatVals : vector of scaling factors (pixelsize)
% epsR : regularizes the division
% doConvolve : Flag that uses a convolved version of f in the nominator, if active.  DEPRECATED! Do not use!
% myReg : Penalty value
% myRegGrad : Gradient
% See also: Verveer et al. Journal of Microscopy, 193, 50-61

function [myReg,myRegGrad]=RegularizeLap27(toRegularize,BetaVals,epsR,doConvolve)
        if nargin < 3
            epsR=2.1;  % If this value is too small the updates can lead to a numerical instability problem
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
            
            aGradLap{1}=(tRL1+tRR1-9*toRegularize)./(BetaVals(1)).^2;
            aGradLap{2}=(tRL2+tRR2-9*toRegularize)./(BetaVals(2)).^2;
            aGradLap{3}=(tRL3+tRR3-9*toRegularize)./(BetaVals(3)).^2;
            drei=(BetaVals(1).^2+BetaVals(2).^2+BetaVals(3).^2);
           
            zweixy=BetaVals(1).^2+BetaVals(2).^2;
            zweixz=BetaVals(1).^2+BetaVals(3).^2;
            zweiyz=BetaVals(2).^2+BetaVals(3).^2;
            aGradLapDU{1}=circshift(tRL1,[0 1 0])/zweixy+circshift(circshift(tRL1,[0 1 0]),[0 0 1])./drei+circshift(circshift(tRL1,[0 1 0]),[0 0 -1])./drei+circshift(tRR1,[0 1 0])/zweixy+circshift(circshift(tRR1,[0 1 0]),[0 0 1])./drei+circshift(circshift(tRR1,[0 1 0]),[0 0 -1])./drei+circshift(circshift(toRegularize,[0 1 0]),[0 0 1])./zweiyz+circshift(circshift(toRegularize,[0 1 0]),[0 0 -1])./zweiyz;
             aGradLapDM{1}=circshift(tRL1,[0 0 1])/zweixz+circshift(tRL1,[0 0 -1])/zweixz+circshift(tRR1,[0 0 1])/zweixz+circshift(tRR1,[0 0 -1])/zweixz;
               aGradLapDO{1}=circshift(tRL1,[0 -1 0])/zweixy+circshift(circshift(tRL1,[0 -1 0]),[0 0 1])/drei+circshift(circshift(tRL1,[0 -1 0]),[0 0 -1])/drei+circshift(tRR1,[0 -1 0])/zweixy+circshift(circshift(tRR1,[0 -1 0]),[0 0 1])/drei+circshift(circshift(tRR1,[0 -1 0]),[0 0 -1])/drei+circshift(circshift(toRegularize,[0 -1 0]),[0 0 1])/zweiyz+circshift(circshift(toRegularize,[0 1 0]),[0 0 -1])/zweiyz;

            %circshift(tRL1,[0 -1 0])+circshift(tRL1,[0 0 -1])+circshift(tRL1,[0 0 1])+circshift(tRR1,[0 1 0])+circshift(tRR1,[0 -1 0])+circshift(tRR1,[0 0 -1])+circshift(tRR1,[0 0 1]);
%             aGradLapD{2}=circshift(tRL2,[1 0 0])+circshift(tRL2,[-1 0 0])+circshift(tRL2,[0 0 -1])+circshift(tRL2,[0 0 1])+circshift(tRR2,[1 0 0])+circshift(tRR2,[-1 0 0])+circshift(tRR2,[0 0 -1])+circshift(tRR2,[0 0 1]);
%             aGradLapD{2}=circshift(tRL3,[1 0 0])+circshift(tRL3,[-1 0 0])+circshift(tRL3,[0 -1 0])+circshift(tRL3,[0 1 0])+circshift(tRR3,[1 0 0])+circshift(tRR3,[-1 0 0])+circshift(tRR3,[0 -1 0])+circshift(tRR3,[0 1 0]);
%             
                     

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
            
            innereAbl=aGradLap{1}.^2+aGradLap{2}.^2+aGradLap{3}.^2+aGradLapDU{1}.^2+aGradLapDM{1}.^2+aGradLapDO{1}.^2;
           
            myReg = sum((innereAbl));
            
           
            myRegGrad1=0;
            myRegGrad2=0;
            myRegGrad3=0;
            myRegGrad4=0;
            myRegGrad5=0;
            myRegGrad6=0;
            myRegGrad7=0;
            
            myRegGrad1 = -18*aGradLap{1}./(BetaVals(1)).^2-18*aGradLap{2}./(BetaVals(2)).^2-18*aGradLap{3}./(BetaVals(3)).^2;
               % circshift(tRL1,[0 1 0])/zweixy+                                circshift(circshift(tRL1,[0 1 0]),[0 0 1])./drei+     circshift(circshift(tRL1,[0 1 0]),[0 0 -1])./drei+       circshift(tRR1,[0 1 0])/zweixy+                    circshift(circshift(tRR1,[0 1 0]),[0 0 1])./drei+      circshift(circshift(tRR1,[0 1 0]),[0 0 -1])./drei+         circshift(circshift(toRegularize,[0 1 0]),[0 0 1])./zweiyz+    circshift(circshift(toRegularize,[0 1 0]),[0 0 -1])./zweiyz;
        
            myRegGrad2=2*circshift(aGradLapDU{1},[0 -1 0])./zweixy+2*circshift(circshift(aGradLapDU{1},[0 -1 0]),[0 0 -1])./drei+2*circshift(circshift(aGradLapDU{1},[0 -1 0]),[0 0 1])./drei+2*circshift(aGradLapDU{1},[0 -1 0])./zweixy+2*circshift(circshift(aGradLapDU{1},[0 -1 0]),[0 0 -1])./drei+circshift(circshift(aGradLapDU{1},[0 -1 0]),[0 0 1])./drei+circshift(circshift(aGradLapDU{1},[0 -1 0]),[0 0 -1])./zweiyz+circshift(circshift(aGradLapDU{1},[0 -1 0]),[0 0 1])./zweiyz;
       %                 circshift(tRL1,[0 0 1])/zweixz+           circshift(tRL1,[0 0 -1])/zweixz+            circshift(tRR1,[0 0 1])/zweixz+             circshift(tRR1,[0 0 -1])/zweixz;
           
            myRegGrad3=2*circshift(aGradLapDM{1},[0 0 -1])./zweixz+2*circshift(aGradLapDM{1},[0 0 1])./zweixz+2*circshift(aGradLapDM{1},[0 0 -1])./zweixz+2*circshift(aGradLapDM{1},[0 0 1])./zweixz;
%             aGradLapDO{1}=circshift(tRL1,[0 -1 0])/zweixy+      circshift(circshift(tRL1,[0 -1 0]),[0 0 1])/drei+            circshift(circshift(tRL1,[0 -1 0]),[0 0 -1])/drei+           circshift(tRR1,[0 -1 0])/zweixy+            circshift(circshift(tRR1,[0 -1 0]),[0 0 1])/drei+            circshift(circshift(tRR1,[0 -1 0]),[0 0 -1])/drei+          circshift(circshift(toRegularize,[0 -1 0]),[0 0 1])/zweiyz+      circshift(circshift(toRegularize,[0 1 0]),[0 0 -1])/zweiyz;

            myRegGrad4=2*circshift(aGradLapDO{1},[0 1 0])./zweixy+2*circshift(circshift(aGradLapDO{1},[0 1 0]),[0 0 -1])./drei+2*circshift(circshift(aGradLapDO{1},[0 1 0]),[0 0 1])./drei+2*circshift(aGradLapDO{1},[0 1 0])./zweixy+2*circshift(circshift(aGradLapDO{1},[0 1 0]),[0 0 -1])./drei+2*circshift(circshift(aGradLapDO{1},[0 1 0]),[0 0 1])./drei+2*circshift(circshift(aGradLapDO{1},[0 1 0]),[0 0 -1])./zweiyz+2*circshift(circshift(aGradLapDO{1},[0 -1 0]),[0 0 1])./zweiyz;
          
            myRegGrad5=2*circshift(aGradLap{1},[1 0 0])./(BetaVals(1)).^2+2*circshift(aGradLap{1},[-1 0 0])./(BetaVals(1)).^2; 
            myRegGrad6=2*circshift(aGradLap{2},[0 1 0])./(BetaVals(2)).^2+2*circshift(aGradLap{2},[0 -1 0])./(BetaVals(2)).^2; 
            myRegGrad7=2*circshift(aGradLap{3},[0 0 1])./(BetaVals(3)).^2+2*circshift(aGradLap{3},[0 0 -1])./(BetaVals(3)).^2; 
                 
                myRegGrad=myRegGrad1+myRegGrad2+myRegGrad3+myRegGrad4+myRegGrad5+myRegGrad6+myRegGrad7;
  end
  %myRegGrad = 2*((aGradL{1}./circshift(nom, [1 0 0]) - aGradR{1}./circshift(nom, [-1 0 0]))/BetaVals(1) + ...
%                         (aGradL{2}./circshift(nom, [0 1 0]) - aGradR{2}./circshift(nom, [0 -1 0]))/BetaVals(2) + ...
%                         (aGradL{3}./circshift(nom, [0 0 1]) - aGradR{3}./circshift(nom, [0 0 -1]))/BetaVals(3) + ...
%                         ((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2) + (aGradL{3} - aGradR{3})/BetaVals(3))./nom) - ...
%                         sign(toRegularizeC).*(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2}) + ...
%                         abssqr(aGradL{3}) + abssqr(aGradR{3}))./nom.^2;