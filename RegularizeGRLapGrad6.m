% [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR,doConvolve) computes Good's roughness regularisation:
% penalty= |Grad(f)|^2 / (|f|+epsR)
% toRegularize : 2D or 3D array to regularize
% BeatVals : vector of scaling factors (pixelsize)
% epsR : regularizes the division
% doConvolve : Flag that uses a convolved version of f in the nominator, if active.  DEPRECATED! Do not use!
% myReg : Penalty value
% myRegGrad : Gradient
% See also: Verveer et al. Journal of Microscopy, 193, 50-61

function [myReg,myRegGrad]=RegularizeGRLapGrad6(toRegularize,BetaVals,epsR,doConvolve)
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
            ab1=BetaVals(1);
            ab2=BetaVals(2);
            ab3=BetaVals(3);
            aGradLap{1}=(tRL1+tRR1-2*toRegularize)./(ab1).^2;
            
            aGradLap{2}=(tRL2+tRR2-2*toRegularize)./(ab2).^2;
            aGradLap{3}=(tRL3+tRR3-2*toRegularize)./(ab3).^2;
            aGradZ{1}=(tRR1 - tRL1)./(2*ab1);
	        aGradZ{2}=(tRR2 - tRL2)./(2*ab2);
            aGradZ{3}=(tRR3 - tRL3)./(2*ab3);
%             aGradZ{1}=(tRR1 - tRL1)/(2*ab1);
%             
% 	    aGradZ{2}=(tRR2 - tRL2)/(2*ab2);
%             aGradZ{3}=(tRR3 - tRL3)/(2*ab3);
%             if doConvolve
%                  toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
%             else
            if doConvolve
                 toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
            else
                 toRegularizeC=toRegularize; 
            end
            epsR=1;
            nom =  abs((aGradZ{1}))+abs((aGradZ{2}))+abs((aGradZ{3}))+epsR;
         %   nom1 = ( abs((aGradZ{1})).*(2*ab2)*(2*ab1)+abs((aGradZ{2})).*(2*ab1)*(2*ab2)+epsR*(2*ab1)*(2*ab2));
            myReg = sum((aGradLap{1}+aGradLap{2}+aGradLap{3})./nom);
% 1d            
%                nom =  abs((aGradZ{1}))+epsR;
%             nom1 = ( abs((aGradZ{1}))*(2*ab1)+epsR*(2*ab1));
%             myReg = sum((aGradLap{1})./nom);
%             
    
            myRegGrad1X=0;
            myRegGrad2X=0;
            myRegGrad3X=0;
            myRegGrad1Y=0;
            myRegGrad2Y=0;
            myRegGrad3Y=0;
            myRegGrad1Z=0;
            myRegGrad2Z=0;
            myRegGrad3Z=0;
            myRegGrad1X2=0;
            myRegGrad1X3=0;
            myRegGrad1Y2=0;
            myRegGrad1Y3=0;
            myRegGrad1Z2=0;
            myRegGrad1Z3=0;
            myRegGrad1X=myRegGrad1X+(-2)./nom./(ab1).^2;
            myRegGrad1X2=myRegGrad1X2+1./(circshift(nom,[-1 0 0]))./(ab1).^2;
            myRegGrad1X3=myRegGrad1X3+1./(circshift(nom,[1 0 0])) ./(ab1).^2;
            
            
            myRegGrad1Y=myRegGrad1Y+(-2)./nom./(ab2).^2;
            myRegGrad1Y2=myRegGrad1Y2+1./(circshift(nom,[0 -1 0]))./(ab2).^2;
            myRegGrad1Y3=myRegGrad1Y3+1./(circshift(nom,[0 1 0])) ./(ab2).^2;
            
            
            
            myRegGrad1Z=myRegGrad1Z+(-2)./nom./(ab3).^2;
            myRegGrad1Z2=myRegGrad1Z2+1./(circshift(nom,[0 0 -1]))./(ab3).^2;
            myRegGrad1Z3=myRegGrad1Z3+1./(circshift(nom,[0 0 1])) ./(ab3).^2;
            
%             
            myRegGrad2X=myRegGrad2X-sign(circshift((aGradZ{1}),[1 0 0])).*(circshift(aGradLap{1}+aGradLap{2}+aGradLap{3},[1 0 0]))./(circshift(nom,[1 0 0]))./((2*ab1));
            
            myRegGrad2Y=myRegGrad2Y-sign(circshift((aGradZ{2}),[0 1 0])).*(circshift(aGradLap{1}+aGradLap{2}+aGradLap{3},[0 1 0]))./(circshift(nom,[0 1 0]))./((2*ab2));
            myRegGrad2Z=myRegGrad2Z-sign(circshift((aGradZ{3}),[0 0 1])).*(circshift(aGradLap{1}+aGradLap{2}+aGradLap{3},[0 0 1]))./(circshift(nom,[0 0 1]))./((2*ab3));

            myRegGrad3X= myRegGrad3X+sign(circshift((aGradZ{1}),[-1 0 0])).*(circshift(aGradLap{1}+aGradLap{2}+aGradLap{3},[-1 0 0]))./(circshift(nom,[-1 0 0]))./((2*ab1));
            
            myRegGrad3Y= myRegGrad3Y+sign(circshift((aGradZ{2}),[0 -1 0])).*(circshift(aGradLap{1}+aGradLap{2}+aGradLap{3},[0 -1 0]))./(circshift(nom,[0 -1 0]))./((2*ab2));
            myRegGrad3Z= myRegGrad3Z+sign(circshift((aGradZ{3}),[0 0 -1])).*(circshift(aGradLap{1}+aGradLap{2}+aGradLap{3},[0 0 -1]))./(circshift(nom,[0 0 -1]))./((2*ab3));

%myRegGrad=myRegGrad2X+myRegGrad2Y+myRegGrad2Z+myRegGrad1X+myRegGrad1Y+myRegGrad1Z+(myRegGrad3X+myRegGrad3Y+myRegGrad3Z)./nom;  
 myRegGrad=myRegGrad1X+myRegGrad1X2+myRegGrad1X3+myRegGrad2X./(circshift(nom,[1 0 0]))+myRegGrad3X./(circshift(nom,[-1 0 0])) ... 
          +myRegGrad1Y+myRegGrad1Y2+myRegGrad1Y3+myRegGrad2Y./(circshift(nom,[0 1 0]))+myRegGrad3Y./(circshift(nom,[0 -1 0])) ...
          +myRegGrad1Z+myRegGrad1Z2+myRegGrad1Z3+myRegGrad2Z./(circshift(nom,[0 0 1]))+myRegGrad3Z./(circshift(nom,[0 0 -1]));     
           
%             myRegGrad = (-2)./nom./(ab1).^2+ (-2)./nom./(ab2).^2+ (-2)./nom./(ab3).^2+1./(circshift(nom{1},[-1 0 0]))./(ab1).^2 ...
%                 +1./(circshift(nom,[1 0 0]))./(ab1).^2+1./(circshift(nom,[0 -1  0]))./(ab2).^2+1./(circshift(nom,[0 1 0]))./(ab2).^2+.1/circshift(nom,[0 0 -1])./(ab3).^2 ...
%                 +1./(circshift(nom,[0 0 1]))./(ab3).^2 ...
%                 -sign(circshift(aGradLap{1},[1 0 0])+circshift(aGradLap{2},[1 0 0])+circshift(aGradLap{3},[1 0  0])).*(circshift(aGradLap{1},[1 0 0]))./(circshift(nom,[1 0  0])).^2 ....
%                 +sign(circshift(aGradLap{1},[-1 0 0])+circshift(aGradLap{2},[-1 0 0])+circshift(aGradLap{3},[-1 0 0])).*(circshift(aGradLap{1},[-1 0 0]))./(circshift(nom,[-1 0 0])).^2 ...
%                 -sign(circshift(aGradLap{1},[0 1 0])+circshift(aGradLap{2},[0 1 0])+circshift(aGradLap{3},[0 1 0])).*(circshift(aGradLap{2},[0 1  0]))./(circshift(nom,[0 1  0])).^2 ....
%                 +sign(circshift(aGradLap{1},[0 -1 0])+circshift(aGradLap{2},[0 -1 0])+circshift(aGradLap{3},[0 -1 0])).*(circshift(aGradLap{2},[0 -1 0]))./(circshift(nom,[0 -1 0])).^2 ...
%                 -sign(circshift(aGradLap{1},[0 0 1])+circshift(aGradLap{2},[0 0 1])+circshift(aGradLap{3},[0 0 1])).*(circshift(aGradLap{3},[0 0 1]))./(circshift(nom,[0 0 1])).^2 ....
%                 +sign(circshift(aGradLap{1},[0 0 -1])+circshift(aGradLap{2},[0 0 -1])+circshift(aGradLap{3},[0 0 -1])).*(circshift(aGradLap{3},[0 0 -1]))./(circshift(nom,[0 0 -1])).^2;



  end
  %myRegGrad = 2*((aGradL{1}./circshift(nom, [1 0 0]) - aGradR{1}./circshift(nom, [-1 0 0]))/ab1 + ...
%                         (aGradL{2}./circshift(nom, [0 1 0]) - aGradR{2}./circshift(nom, [0 -1 0]))/ab2 + ...
%                         (aGradL{3}./circshift(nom, [0 0 1]) - aGradR{3}./circshift(nom, [0 0 -1]))/ab3 + ...
%                         ((aGradL{1} - aGradR{1})/ab1 + (aGradL{2} - aGradR{2})/ab2 + (aGradL{3} - aGradR{3})/ab3)./nom) - ...
%                         sign(toRegularizeC).*(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2}) + ...
%                         abssqr(aGradL{3}) + abssqr(aGradR{3}))./nom.^2;