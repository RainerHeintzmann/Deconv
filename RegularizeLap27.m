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
            epsR=0.1;  % If this value is too small the updates can lead to a numerical instability problem
        end
        if nargin < 4
            doConvolve=0;   % If Active: This Trick devides by a better estimate of the local intensity.
        end
        if doConvolve~=0
            error('doConvolve option is deprecated. Do not use.');
        end
       
        
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

            
            innereAbl=aGradLap{1}+aGradLap{2}+aGradLap{3}+aGradLapDU{1}+aGradLapDM{1}+aGradLapDO{1};
           
            myReg = sum((innereAbl).^2);
            
           
            myRegGrad1=0;
            myRegGrad2=0;
            myRegGrad3=0;
            myRegGrad4=0;
            myRegGrad5=0;
            myRegGrad6=0;
            myRegGrad7=0;

            myRegGrad1 = -18*innereAbl./(BetaVals(1)).^2-18*innereAbl./(BetaVals(2)).^2-18*innereAbl./(BetaVals(3)).^2;
                           % circshift(tRL1,[0 1 0])/zweixy+                            circshift(circshift(tRL1,[0 1 0]),[0 0 1])./drei+                               circshift(circshift(tRL1,[0 1 0]),[0 0 -1])./drei+                        circshift(tRR1,[0 1 0])/zweixy+                               circshift(circshift(tRR1,[0 1 0]),[0 0 1])./drei+                        circshift(circshift(tRR1,[0 1 0]),[0 0 -1])./drei+                         circshift(circshift(toRegularize,[0 1 0]),[0 0 1])./zweiyz+    circshift(circshift(toRegularize,[0 1 0]),[0 0 -1])./zweiyz;
              myRegGrad2=2*circshift(circshift(innereAbl,[0 -1 0]),[-1 0 0])./zweixy+2*circshift(circshift(circshift(innereAbl,[0 -1 0]),[0 0 -1]),[-1 0 0])./drei+2*circshift(circshift(circshift(innereAbl,[0 -1 0]),[0 0 1]),[-1 0 0])./drei+2*circshift(circshift(innereAbl,[0 -1 0]),[1 0 0])./zweixy+2*circshift(circshift(circshift(innereAbl,[0 -1 0]),[0 0 -1]),[1 0 0])./drei+circshift(circshift(circshift(innereAbl,[0 -1 0]),[0 0 1]),[1 0 0])./drei+circshift(circshift(innereAbl,[0 -1 0]),[0 0 -1])./zweiyz+circshift(circshift(innereAbl,[0 -1 0]),[0 0 1])./zweiyz;
           %                 circshift(tRL1,[0 0 1])/zweixz+                          circshift(tRL1,[0 0 -1])/zweixz+                              circshift(tRR1,[0 0 1])/zweixz+                            circshift(tRR1,[0 0 -1])/zweixz;
           
            myRegGrad3=2*circshift(circshift(innereAbl,[0 0 -1]),[-1 0 0])./zweixz+2*circshift(circshift(innereAbl,[0 0 1]),[-1 0 0])./zweixz+2*circshift(circshift(innereAbl,[0 0 -1]),[1 0 0])./zweixz+2*circshift(circshift(innereAbl,[0 0 1]),[1 0 0])./zweixz;
%             aGradLapDO{1}=circshift(tRL1,[0 -1 0])/zweixy+                  circshift(circshift(tRL1,[0 -1 0]),[0 0 1])/drei+                          circshift(circshift(tRL1,[0 -1 0]),[0 0 -1])/drei+                             circshift(tRR1,[0 -1 0])/zweixy+                              circshift(circshift(tRR1,[0 -1 0]),[0 0 1])/drei+                                circshift(circshift(tRR1,[0 -1 0]),[0 0 -1])/drei+               circshift(circshift(toRegularize,[0 -1 0]),[0 0 1])/zweiyz+circshift(circshift(toRegularize,[0 1 0]),[0 0 -1])/zweiyz;
         myRegGrad4=2*circshift(circshift(innereAbl,[0 1 0]),[-1 0 0])./zweixy+2*circshift(circshift(circshift(innereAbl,[0 1 0]),[0 0 -1]),[-1 0 0])./drei+2*circshift(circshift(circshift(innereAbl,[0 1 0]),[0 0 1]),[-1 0 0])./drei+2*circshift(circshift(innereAbl,[0 1 0]),[1 0 0])./zweixy+2*circshift(circshift(circshift(innereAbl,[0 1 0]),[0 0 -1]),[1 0 0])./drei+2*circshift(circshift(circshift(innereAbl,[0 1 0]),[0 0 1]),[1 0 0])./drei+2*circshift(circshift(innereAbl,[0 1 0]),[0 0 -1])./zweiyz+2*circshift(circshift(innereAbl,[0 -1 0]),[0 0 1])./zweiyz;
          
            myRegGrad5=2*circshift(innereAbl,[1 0 0])./(BetaVals(1)).^2+2*circshift(innereAbl,[-1 0 0])./(BetaVals(1)).^2; 
            myRegGrad6=2*circshift(innereAbl,[0 1 0])./(BetaVals(2)).^2+2*circshift(innereAbl,[0 -1 0])./(BetaVals(2)).^2; 
            myRegGrad7=2*circshift(innereAbl,[0 0 1])./(BetaVals(3)).^2+2*circshift(innereAbl,[0 0 -1])./(BetaVals(3)).^2; 
                 
                myRegGrad=myRegGrad1+myRegGrad2+myRegGrad3+myRegGrad4+myRegGrad5+myRegGrad6+myRegGrad7;
  end
               