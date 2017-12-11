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
        nom = epsR + abs(toRegularize);
        
        myRegGrad_AbsSqr = 0;  myRegGrad = 0; myRegGrad_B = 0;
        numdims=ndims(toRegularize);
        for d=1:numdims     % This algorithm is n-dimensional. If a particular dimension should not be regularized: choose BetaVals(d)=0 or NaN
            if ~isnan(BetaVals(d)) && BetaVals(d)~=0
                myrotshift = zeros(1,numdims); myrotshift(d) = 1.0;
                
                tRL=circshift(toRegularize,myrotshift);             % cyclic rotation
                aGradL=(toRegularize - tRL)/BetaVals(d);      % left-sided gradient X
                clear tRL;
                tRR=circshift(toRegularize,-myrotshift);
                aGradR=(tRR - toRegularize)/BetaVals(d);      % right-sided gradient X
                clear tRR;
                myRegGrad_AbsSqr = myRegGrad_AbsSqr  + abssqr(aGradL) + abssqr(aGradR);  % This is the sum that needs to be devided by the common denominator later
                myRegGrad = myRegGrad  + (aGradL./circshift(nom, myrotshift) - aGradR./circshift(nom, -myrotshift))/BetaVals(d);
                myRegGrad_B = myRegGrad_B + (aGradL - aGradR)/BetaVals(d);  % This is the sum that needs to be devided by the common denominator later
            end
        end
        clear aGradL; clear aGradR;        
        myRegGrad_AbsSqr  = myRegGrad_AbsSqr ./ nom;
        myReg = sum(myRegGrad_AbsSqr);   % is already devided by "nom"
        myRegGrad = 2*(myRegGrad + myRegGrad_B./nom);
        clear myRegGrad_B;
        myRegGrad = myRegGrad - sign(toRegularize).*myRegGrad_AbsSqr ./ nom;  % myRegGrad_AbsSqr is already once devided by "nom", now a second time
end