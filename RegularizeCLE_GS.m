% [myReg,myRegGrad]=RegularizeGS(toRegularize,BetaVals) computes a Gradient squared regularisation. This is similar to Good's roughness but not the same
% penalty=|grad(f)|^2
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient
% g=gradient(EM);
% RefImgX=g{1}.^2;  % close to zero means high regularization, high value mean no regularization.
% RefImgY=g{2}.^2;  % close to zero means high regularization, high value mean no regularization.
% RefImgX=1-RefImgX/max(RefImgX);
% RefImgY=1-RefImgY/max(RefImgY);


function [myReg,myRegGrad]=RegularizeCLE_GS(toRegularize,BetaVals,RefImgX,RefImgY)   
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
                myReg = sum((abssqr(aGradL{1})+abssqr(aGradR{1}))./RefImgX + (abssqr(aGradL{2})+abssqr(aGradR{2}))./RefImgY);
            % gradient is tested
                myRegGrad = 4*(((aGradL{1} - aGradR{1})/BetaVals(1))./RefImgX + ((aGradL{2} - aGradR{2})/BetaVals(2))./RefImgY);
                  
        elseif ndims(toRegularize) == 3
            aGradL{1}=(toRegularize - circshift(toRegularize,[1 0 0]))/BetaVals(1);  % cyclic rotation
            aGradL{2}=(toRegularize - circshift(toRegularize,[0 1 0]))/BetaVals(2);  % cyclic rotation
            aGradL{3}=(toRegularize - circshift(toRegularize,[0 0 1]))/BetaVals(3);  % cyclic rotation
            aGradR{1}=(circshift(toRegularize,[-1 0 0]) - toRegularize)/BetaVals(1);
            aGradR{2}=(circshift(toRegularize,[0 -1 0]) - toRegularize)/BetaVals(2);
            aGradR{3}=(circshift(toRegularize,[0 0 -1]) - toRegularize)/BetaVals(3);
            
            myReg = sum((abssqr(aGradL{1}))./RefImgX + (abssqr(aGradL{2}))./RefImgY + abssqr(aGradL{3}) + (abssqr(aGradR{1}))./RefImgX + (abssqr(aGradR{2}))./RefImgY + abssqr(aGradR{3}));
            myRegGrad = 4*(((aGradL{1} - aGradR{1})/BetaVals(1))./RefImgX + ((aGradL{2} - aGradR{2})/BetaVals(2))./RefImgY + (aGradL{3} - aGradR{3})/BetaVals(3));       
        else % 1-D
            aGradL=(toRegularize - circshift(toRegularize,1))/BetaVals(1);  
            aGradR=(circshift(toRegularize,-1) - toRegularize)/BetaVals(1);  % cyclic rotation
            myReg = sum((abssqr(aGradL) + abssqr(aGradR))./RefImgX);
            myRegGrad = 2*(((aGradL - aGradR)/BetaVals(1) + (aGradL - aGradR)/BetaVals(1))./RefImgX);
        end            
        