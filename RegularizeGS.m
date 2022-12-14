% [myReg,myRegGrad]=RegularizeGS(toRegularize,BetaVals) computes a Gradient squared regularisation. This is similar to Good's roughness but not the same
% penalty=|grad(f)|^2
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient

function [myReg,myRegGrad]=RegularizeGS(toRegularize,BetaVals)
        % aGrad=gradient(aRecon)
        % myHessian = hessian(aRecon);
        % myRegGrad = -2*real(myHessian{1,1}+myHessian{2,2}+2*myHessian{1,2}); 
        if (ndims(toRegularize) == 2) || (size(toRegularize,3) == 1)
            if (ndims(toRegularize) == 2)
                aGradL{1}=(toRegularize - circshift(toRegularize,[1 0]))/BetaVals(1);  % cyclic rotation
                aGradL{2}=(toRegularize - circshift(toRegularize,[0 1]))/BetaVals(2);  % cyclic rotation
                aGradR{1}=(circshift(toRegularize,[-1 0]) - toRegularize)/BetaVals(1);
                aGradR{2}=(circshift(toRegularize,[0 -1]) - toRegularize)/BetaVals(2);
            else (size(toRegularize,3) == 1)
                aGradL{1}=(toRegularize - circshift(toRegularize,[1 0 0]))/BetaVals(1);  % cyclic rotation
                aGradL{2}=(toRegularize - circshift(toRegularize,[0 1 0]))/BetaVals(2);  % cyclic rotation
                aGradR{1}=(circshift(toRegularize,[-1 0 0]) - toRegularize)/BetaVals(1);
                aGradR{2}=(circshift(toRegularize,[0 -1 0]) - toRegularize)/BetaVals(2);
            end
            myReg = sum(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradR{1}) + abssqr(aGradR{2}));
            % gradient is tested
%             myRegGrad = 2*((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2) + ...
%                         (aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2));
            myRegGrad = 4*((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2));
        elseif ndims(toRegularize) == 3
            aGradL{1}=(toRegularize - circshift(toRegularize,[1 0 0]))/BetaVals(1);  % cyclic rotation
            aGradL{2}=(toRegularize - circshift(toRegularize,[0 1 0]))/BetaVals(2);  % cyclic rotation
            aGradL{3}=(toRegularize - circshift(toRegularize,[0 0 1]))/BetaVals(3);  % cyclic rotation
            aGradR{1}=(circshift(toRegularize,[-1 0 0]) - toRegularize)/BetaVals(1);
            aGradR{2}=(circshift(toRegularize,[0 -1 0]) - toRegularize)/BetaVals(2);
            aGradR{3}=(circshift(toRegularize,[0 0 -1]) - toRegularize)/BetaVals(3);
            myReg = sum(abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradL{3}) + abssqr(aGradR{1}) + abssqr(aGradR{2}) + abssqr(aGradR{3}));
%             myRegGrad = 2*((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2) + (aGradL{3} - aGradR{3})/BetaVals(3) + ...
%                         (aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2) + (aGradL{3} - aGradR{3})/BetaVals(3));
            myRegGrad = 4*((aGradL{1} - aGradR{1})/BetaVals(1) + (aGradL{2} - aGradR{2})/BetaVals(2) + (aGradL{3} - aGradR{3})/BetaVals(3));
        else % 1-D
            aGradL=(toRegularize - circshift(toRegularize,1))/BetaVals(1);  
            aGradR=(circshift(toRegularize,-1) - toRegularize)/BetaVals(1);  % cyclic rotation
            myReg = sum(abssqr(aGradL) + abssqr(aGradR));
            % myRegGrad = 2*((aGradL - aGradR)/BetaVals(1) + (aGradL - aGradR)/BetaVals(1));
            myRegGrad = 4*((aGradL - aGradR)/BetaVals(1));
        end
        