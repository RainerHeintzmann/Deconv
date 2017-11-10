% [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR,doConvolve) computes Good's roughness regularisation:
% penalty= |Grad(f)|^2 / (|f|+epsR)
% toRegularize : 2D or 3D array to regularize
% BeatVals : vector of scaling factors (pixelsize)
% epsR : regularizes the division
% doConvolve : Flag that uses a convolved version of f in the nominator, if active.  DEPRECATED! Do not use!
% myReg : Penalty value
% myRegGrad : Gradient
% See also: Verveer et al. Journal of Microscopy, 193, 50-61

function [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR,doConvolve)
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
            clear tRR3 tRR2 tRR1 tRL1 tRL2 tRL3
%             aGradZ{1}=(tRR1 - tRL1)/(2*ab1);
%             
% 	    aGradZ{2}=(tRR2 - tRL2)/(2*ab2);
%             aGradZ{3}=(tRR3 - tRL3)/(2*ab3);
%             if doConvolve
%                  toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
%             else
%            if doConvolve
%                toRegularizeC=(toRegularize+tRL1+tRR1+tRL2+tRR2+tRL3+tRR3)/7;  % blurs the nominator
%           else
%                 toRegularizeC=toRegularize; 
%           end
            
           above= (aGradLap{1}+aGradLap{2}+aGradLap{3});
		   clear aGradLap
          
            nom =  ((aGradZ{1})).^2+((aGradZ{2})).^2+((aGradZ{3})).^2+(toRegularize).^2+epsR;
         %   nom1 = ( abs((aGradZ{1})).*(2*ab2)*(2*ab1)+abs((aGradZ{2})).*(2*ab1)*(2*ab2)+epsR*(2*ab1)*(2*ab2));
            myReg = sum((above).^2./nom);
        
        myRegGrad1Above2 = 0;  
        numdims=ndims(toRegularize);
		myRegGrad1Above=0;
        for d=1:numdims     % This algorithm is n-dimensional. If a particular dimension should not be regularized: choose BetaVals(d)=0 or NaN
            if ~isnan(BetaVals(d)) && BetaVals(d)~=0
                myrotshift = zeros(1,numdims); myrotshift(d) = 1.0;
                
                
				myRegGrad1Above=myRegGrad1Above+2*(above).*(-2)./nom./(BetaVals(d)).^2;
				myRegGrad1Above=myRegGrad1Above+2*(circshift(above,myrotshift))./(circshift(nom,myrotshift)) ./(BetaVals(d)).^2;
                     myrotshift(d) = -1.0;  
		   myRegGrad1Above=myRegGrad1Above+2*(circshift(above,myrotshift))./(circshift(nom,myrotshift))./(BetaVals(d)).^2;
            
  
            end
        end
		for d=1:numdims     % This algorithm is n-dimensional. If a particular dimension should not be regularized: choose BetaVals(d)=0 or NaN
            if ~isnan(BetaVals(d)) && BetaVals(d)~=0
                myrotshift = zeros(1,numdims); myrotshift(d) = 1.0;
                
				 
  myRegGrad1Above2=myRegGrad1Above2-2*(circshift((aGradZ{d}),myrotshift)).*(circshift(above,myrotshift)).^2./(circshift(nom,myrotshift)).^2./((2*BetaVals(d)));
                                myrotshift(d) = -1.0;  
		 myRegGrad1Above2= myRegGrad1Above2+2*(circshift((aGradZ{d}),myrotshift)).*(circshift(above,myrotshift)).^2./(circshift(nom,myrotshift)).^2./((2*BetaVals(d)));
               
  
            end
        end
		 myRegGrad1Above2=  myRegGrad1Above2-2*(toRegularize).*(above).^2./(nom)./nom;
       
        myRegGrad =myRegGrad1Above+ myRegGrad1Above2;  % myRegGrad_AbsSqr is already once devided by "nom", now a second time
