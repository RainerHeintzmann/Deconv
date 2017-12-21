% [myReg,myRegGrad]=RegularizeLap(toRegularize,BetaVals,epsR,doConvolve) computes regularisation:
% penalty= (Laplace(f))^2 
% toRegularize : 2D or 3D array to regularize
% BeatVals : vector of scaling factors (pixelsize)
% epsR : regularizes the division
% doConvolve : Flag that uses a convolved version of f in the nominator, if active.  DEPRECATED! Do not use!
% myReg : Penalty value
% myRegGrad : Gradient
% Written by Felix Mucha: f.mucha@web.de

function [myReg,myRegGrad]=RegularizeLap(toRegularize,BetaVals,doConvolve)
       
        if nargin < 3
            doConvolve=0;   % If Active: This Trick devides by a better estimate of the local intensity.
        end
        if doConvolve~=0
            error('doConvolve option is deprecated. Do not use.');
        end
       
        
            tRL1=circshift(toRegularize,[1 0 0]);tRR1=circshift(toRegularize,[-1 0 0]);
            tRL2=circshift(toRegularize,[0 1 0]);tRR2=circshift(toRegularize,[0 -1 0]);
            tRL3=circshift(toRegularize,[0 0 1]);tRR3=circshift(toRegularize,[0 0 -1]);
            aGradLap{1}=(tRL1+tRR1-2*toRegularize)/(BetaVals(1)).^2;
            aGradLap{2}=(tRL2+tRR2-2*toRegularize)/(BetaVals(2)).^2;
            aGradLap{3}=(tRL3+tRR3-2*toRegularize)/(BetaVals(3)).^2;
            clear tRL1 tRL2 tRL3 tRR1 tRR2 tRR3

            myReg = sum((aGradLap{1}+aGradLap{2}+aGradLap{3}).^2);
            above= (aGradLap{1}+aGradLap{2}+aGradLap{3});
            clear aGradLap
            numdims=ndims(toRegularize);
            clear toRegularize
		    myRegGrad=0;
       for d=1:numdims     % This algorithm is n-dimensional. If a particular dimension should not be regularized: choose BetaVals(d)=0 or NaN
            if ~isnan(BetaVals(d)) && BetaVals(d)~=0
                myrotshift = zeros(1,numdims); myrotshift(d) = 1.0;
                
                
				myRegGrad=myRegGrad+2*(above).*(-2)/(BetaVals(d)).^2;
				myRegGrad=myRegGrad+2*(circshift(above,myrotshift)) ./(BetaVals(d)).^2;
                myrotshift(d) = -1.0;  
		        myRegGrad=myRegGrad+2*(circshift(above,myrotshift))./(BetaVals(d)).^2;
            end
       end            
  end
  