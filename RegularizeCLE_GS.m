% [myReg,myRegGrad]=RegularizeCLE_GS(toRegularize,BetaVals) computes a guided Gradient squared regularisation.
% penalty=|grad(f)|^2/RefImg
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient
% RefL and RefR should be calculated in the following way:
% RefL{1}= EM - circshift(EM,[1,0]);
% RefL{2}= EM - circshift(EM,[0,1]);
% RefR{1}= circshift(EM,[-1,0] - EM);
% RefR{2}= circshift(EM,[-0,1] - EM);
% RefL=alpha*(|RefL{1}|^2+|RefL{2}|^2)+gamma
% RefR=alpha*(|RefR{1}|^2+|RefR{2}|^2)+gamma

% RefImgX=RefL;  % close to zero means high regularization, high value mean no regularization.
% RefImgY=RefR;  % close to zero means high regularization, high value mean no regularization.


function [myReg,myRegGrad]=RegularizeCLE_GS(toRegularize,BetaVals,RefImgX,RefImgY)   
     %both RefImgX and RefImgY are EM image
     %the gradient of EM is calculated by circshift(idealy not calculated
     %during the iterations)
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
%                 myReg = sum((abssqr(aGradL{1}) + abssqr(aGradL{2}))./RefImgX + (abssqr(aGradR{1})  + abssqr(aGradR{2}))./RefImgY);
                myReg = sum((abssqr(aGradL{1}) + abssqr(aGradL{2}))./RefImgX + (abssqr(aGradR{1})  + abssqr(aGradR{2}))./RefImgY );
            % gradient is tested
                myRegGrad = 4*((aGradL{1}/BetaVals(1) + aGradL{2}/BetaVals(2))./RefImgX - (aGradR{1}/BetaVals(1)+ aGradR{2}/BetaVals(2))./RefImgY - ...
                    aGradR{1}/BetaVals(1)./circshift(RefImgX,[-1 0]) - aGradR{2}/BetaVals(2)./circshift(RefImgX, [0,-1]) + ...
                    aGradL{1}/BetaVals(1)./circshift(RefImgY,[1,0]) + aGradL{2}/BetaVals(2)./circshift(RefImgY,[0,1]));
                  
        elseif ndims(toRegularize) == 3
            aGradL{1}=(toRegularize - circshift(toRegularize,[1 0 0]))/BetaVals(1);  % cyclic rotation
            aGradL{2}=(toRegularize - circshift(toRegularize,[0 1 0]))/BetaVals(2);  % cyclic rotation
            aGradL{3}=(toRegularize - circshift(toRegularize,[0 0 1]))/BetaVals(3);  % cyclic rotation
            aGradR{1}=(circshift(toRegularize,[-1 0 0]) - toRegularize)/BetaVals(1);
            aGradR{2}=(circshift(toRegularize,[0 -1 0]) - toRegularize)/BetaVals(2);
            aGradR{3}=(circshift(toRegularize,[0 0 -1]) - toRegularize)/BetaVals(3);
            
            myReg = sum((abssqr(aGradL{1}) + abssqr(aGradL{2}) + abssqr(aGradL{3}))./RefImgX + (abssqr(aGradR{1}) + abssqr(aGradR{2}) + abssqr(aGradR{3}))./RefImgY);
            myRegGrad = 4*((aGradL{1}/BetaVals(1) + aGradL{2}/BetaVals(2) + aGradL{3}/BetaVals(3))./RefImgX -...
                           (aGradR{1}/BetaVals(1) + aGradR{2}/BetaVals(2) + aGradR{3}/BetaVals(3))./RefImgY - ...
                            aGradR{1}/BetaVals(1)./circshift(RefImgX,[-1,0,0]) - aGradR{2}/BetaVals(2)./circshift(RefImgX, [0,-1,0]) - ...
                            aGradR{3}/BetaVals(3)./circshift(RefImgX, [0,0,-1]) + aGradL{1}/BetaVals(1)./circshift(RefImgY, [1,0,0]) + ...
                            aGradL{2}/BetaVals(2)./circshift(RefImgY,[0,1,0]) + aGradL{3}/BetaVals(3)./circshift(RefImgY, [0,0,1]));       
        else % 1-D
            aGradL=(toRegularize - circshift(toRegularize,1))/BetaVals(1);  
            aGradR=(circshift(toRegularize,-1) - toRegularize)/BetaVals(1);  % cyclic rotation
            myReg = sum(abssqr(aGradL)./RefImgX + abssqr(aGradR)./RefImgY);
            myRegGrad = 4*((aGradL/BetaVals(1))./RefImgX - (aGradR/BetaVals(1))./RefImgY - ...
                aGradR/BetaVals(1)./circshift(RefImgX, -1) + aGradL/BetaVals(1)./circshift(RefImgY,1));
        end            
        