% [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals) computes Good's roughness regularisation:
% penalty= |Grad(f)|^2 / (|f|+epsR)
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient

function [myReg,myRegGrad]=RegularizeGR(toRegularize,BetaVals,epsR)
        % aGrad=gradient(aRecon)
        % myHessian = hessian(aRecon);
        % myRegGrad = -2*real(myHessian{1,1}+myHessian{2,2}+2*myHessian{1,2}); 
        if nargin < 3
            epsR=1e-4;
        end
        if (ndims(toRegularize) == 2) || (size(toRegularize,3) == 1)
            if (ndims(toRegularize) == 2)
                aGrad{1}=(circshift(toRegularize,[1 0])-circshift(toRegularize,[-1 0]))/(2*BetaVals(1));  % cyclic rotation
                aGrad{2}=(circshift(toRegularize,[0 1])-circshift(toRegularize,[0 -1]))/(2*BetaVals(2));  % cyclic rotation
                myRegGrad = ((1/(2*BetaVals(1)*BetaVals(1)))*(2*toRegularize-circshift(toRegularize,[-2 0]) - circshift(toRegularize,[2 0])) + ...
                    (1/(2*BetaVals(2)*BetaVals(2)))*(2*toRegularize-circshift(toRegularize,[0 -2]) - circshift(toRegularize,[0 2]))) ./ (epsR+abs(toRegularize)) - ...
                    sign(toRegularize)*((aGrad{1}).^2+(aGrad{2}).^2)./((epsR+abs(toRegularize)).^2);
            else (size(toRegularize,3) == 1)
                aGrad{1}=(circshift(toRegularize,[1 0 0])-circshift(toRegularize,[-1 0 0]))/2;  % cyclic rotation
                aGrad{2}=(circshift(toRegularize,[0 1 0])-circshift(toRegularize,[0 -1 0]))/2;  % cyclic rotation
                myRegGrad = ((1/(2*BetaVals(1)*BetaVals(1)))*(2*toRegularize-circshift(toRegularize,[-2 0 0]) - circshift(toRegularize,[2 0 0])) + ...
                    (1/(2*BetaVals(2)*BetaVals(2)))*(2*toRegularize-circshift(toRegularize,[0 -2 0]) - circshift(toRegularize,[0 2 0])))  ./ (epsR+abs(toRegularize)) - ...
                    sign(toRegularize)*((aGrad{1}).^2+(aGrad{2}).^2)./((epsR+abs(toRegularize)).^2);
            end
        
            % myReg = sum(BetaVals(1)*aGrad{1} .* aGrad{1} + BetaVals(2)*aGrad{2} .* aGrad{2});
            myReg = sum((abssqr(aGrad{1}) + abssqr(aGrad{2})) ./ (epsR+abs(toRegularize)));
        elseif ndims(toRegularize) == 3
            aGrad{1}=(circshift(toRegularize,[1 0 0])-circshift(toRegularize,[-1 0 0]))/(2*BetaVals(1));  % cyclic rotation
            aGrad{2}=(circshift(toRegularize,[0 1 0])-circshift(toRegularize,[0 -1 0]))/(2*BetaVals(2));  % cyclic rotation
            aGrad{3}=(circshift(toRegularize,[0 0 1])-circshift(toRegularize,[0 0 -1]))/(2*BetaVals(3));  % cyclic rotation
            myReg = sum((abssqr(aGrad{1}) + abssqr(aGrad{2}) + abssqr(aGrad{3})) ./ (epsR+abs(toRegularize)));
            myRegGrad = ((1/(2*BetaVals(1)*BetaVals(1)))*(2*toRegularize-circshift(toRegularize,[-2 0 0])-circshift(toRegularize,[2 0 0])) + ...
                        (1/(2*BetaVals(2)*BetaVals(2)))*(2*toRegularize-circshift(toRegularize,[0 -2 0])-circshift(toRegularize,[0 2 0])) + ...
                        (1/(2*BetaVals(3)*BetaVals(3)))*(2*toRegularize-circshift(toRegularize,[0 0 -2])-circshift(toRegularize,[0 0 2]))) ./ (epsR+abs(toRegularize)) - ...
                        sign(toRegularize)*((aGrad{1}).^2+(aGrad{2}).^2+(aGrad{3}).^2)./((epsR+abs(toRegularize)).^2);
        else % 1-D
            aGrad{1}=(circshift(toRegularize,1)-circshift(toRegularize,-1))/(2*BetaVals(1));  % cyclic rotation
            myReg = sum(abssqr(aGrad{1}) ./ (epsR+abs(toRegularize)));
            myRegGrad = ((1/(2*BetaVals(1)*BetaVals(1)))*(2*toRegularize-circshift(toRegularize,-2) - circshift(toRegularize,2)))  ./ (epsR+abs(toRegularize)) - ...
                        sign(toRegularize)*((aGrad{1}).^2)./((epsR+abs(toRegularize)).^2);
        end 
