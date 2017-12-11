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
      
       
        
            tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
            tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
            tRL3=circshift(toRegularize,[0 0 1]);tRR3=circshift(toRegularize,[0 0 -1]);
            
    	    aGrad{1}=(tRR1 - tRL1)./(2*BetaVals(1));
            aGrad{2}=(tRR2 - tRL2)./(2*BetaVals(2));
            aGrad{3}=(tRR3 - tRL3)./(2*BetaVals(3));
     
            clear tRR1 tRR2 tRR3 tRL1 tRL2 tRL3
            nom = epsR + abs(toRegularize);
            myReg = sum((((abssqr( aGrad{1})+abssqr( aGrad{2})+abssqr( aGrad{3}))./nom))  );
            above=(abssqr( aGrad{1})+abssqr( aGrad{2})+abssqr( aGrad{3}));
                       
          
          
           numdims=ndims(toRegularize);
           clear toRegularize
               
		myRegGrad=0;
       for d=1:numdims     % This algorithm is n-dimensional. If a particular dimension should not be regularized: choose BetaVals(d)=0 or NaN
            if ~isnan(BetaVals(d)) && BetaVals(d)~=0
                myrotshift = zeros(1,numdims); myrotshift(d) = 1.0;
                
                
				myRegGrad=myRegGrad+2*(((circshift( aGrad{d}, myrotshift))))./(2*BetaVals(d))./(circshift(nom, myrotshift));
				
                     myrotshift(d) = -1.0;  
		  myRegGrad=myRegGrad-2*(((circshift( aGrad{d}, myrotshift))))./(2*BetaVals(d))./(circshift(nom, myrotshift));
            
  
            end
       end
        clear aGrad
        myRegGrad=myRegGrad-sign(nom).* (above)./(nom)./nom;
            
            
            
         
end
 