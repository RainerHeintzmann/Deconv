% [myReg,myRegGrad]=RegularizeIG(toRegularize,BetaVals) computes the
% intensity-guided deconvolution for CLEM
% penalty=f/EM
% toRegularize : 2D or 3D array to regularize
% myReg : Penalty value
% myRegGrad : Gradient (a constant value)
% g=gradient(EM);
% RefImg=EM ;  % close to zero means high regularization, high value mean no regularization.



function [myReg,myRegGrad]=RegularizeIG(toRegularize,BetaVals,RefImg,ep)  
             RefImg = RefImg + ep;% fengjiao 25.05.2020
             myReg = sum(toRegularize./RefImg);
             myRegGrad = 1/RefImg;                
end            
        