% [myReg,myRegGrad]=RegularizeGRLapDivReg(toRegularize,BetaVals,epsR,doConvolve) computes regularisation:
% penalty= |Lap(f)|^2 / (|grad(f+v)|+|f|+epsR)
% toRegularize : 2D or 3D array to regularize
% BeatVals : vector of scaling factors (pixelsize)
% epsR,v : regularizes the division
% doConvolve : Flag that uses a convolved version of f in the nominator, if active.  DEPRECATED! Do not use!
% myReg : Penalty value
% myRegGrad : Gradient
% Written by Felix Mucha: f.mucha@web.de

function [myReg,myRegGrad]=RegularizeGRLapDivReg(toRegularize,BetaVals,epsR,doConvolve)
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

         
           above= (aGradLap{1}+aGradLap{2}+aGradLap{3});
		   clear aGradLap
           v=0.001;

            nom =  (abs(aGradZ{1}+v)+abs(aGradZ{2}+v)+abs(aGradZ{3}+v)+abs(toRegularize)).^2+epsR;
            nomo =  ((abs(aGradZ{1}+v))+(abs(aGradZ{2}+v))+(abs(aGradZ{3}+v))+abs(toRegularize));
        
            myReg = sum((above).^2./nom);
        
        myRegGrad = 0;  
        numdims=ndims(toRegularize);

        for d=1:numdims     % This algorithm is n-dimensional. If a particular dimension should not be regularized: choose BetaVals(d)=0 or NaN
            if ~isnan(BetaVals(d)) && BetaVals(d)~=0
                myrotshift = zeros(1,numdims); myrotshift(d) = 1.0;
                
                
                myRegGrad=myRegGrad-2.*sign(circshift((aGradZ{d}+v),myrotshift)).*(circshift((nomo),myrotshift)).*(circshift(above,myrotshift)).^2./(circshift(nom,myrotshift)).^2./((2*BetaVals(d)));
                myrotshift(d) = -1.0;
                myRegGrad= myRegGrad+2.*sign(circshift((aGradZ{d}+v),myrotshift)).*(circshift((nomo),myrotshift)).*(circshift(above,myrotshift)).^2./(circshift(nom,myrotshift)).^2./((2*BetaVals(d)));
                
            end
        end
        clear aGradZ
        for d=1:numdims     % This algorithm is n-dimensional. If a particular dimension should not be regularized: choose BetaVals(d)=0 or NaN
            if ~isnan(BetaVals(d)) && BetaVals(d)~=0
                myrotshift = zeros(1,numdims); myrotshift(d) = 1.0;
                
                
				myRegGrad=myRegGrad+2*(above).*(-2)./nom./(BetaVals(d)).^2;
				myRegGrad=myRegGrad+2*(circshift(above,myrotshift))./(circshift(nom,myrotshift)) ./(BetaVals(d)).^2;
                     myrotshift(d) = -1.0;  
		   myRegGrad=myRegGrad+2*(circshift(above,myrotshift))./(circshift(nom,myrotshift))./(BetaVals(d)).^2;
            
  
            end
        end

		 myRegGrad=  myRegGrad-2.*sign(toRegularize).*(((nomo))).*(above).^2./(nom)./nom;
end