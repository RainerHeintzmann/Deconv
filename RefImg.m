function [RefImgX,RefImgY]=RefImg(EM,x,BetaVals)
if nargin < 4
        BetaVals=[1 1 1];
end 
if nargin < 2
        x=2;
end 
%% guide image
if (ndims(EM) == 2) || (size(EM,3) == 1)% EM - 2D
        if (size(EM,3) == 1)
            EM=squeeze(EM);
        end
        RefL{1}=(EM-circshift(EM,[1,0]))/BetaVals(1);
        RefL{2}=(EM-circshift(EM,[0,1]))/BetaVals(2);
        RefR{1}=(circshift(EM,[-1,0])-EM)/BetaVals(1);
        RefR{2}=(circshift(EM,[0,-1])-EM)/BetaVals(2);
        RefImgX=sqrt(abssqr(RefL{1})+abssqr(RefL{2}));
        RefImgY=sqrt(abssqr(RefR{1})+abssqr(RefR{2}));
        RefImgX=(RefImgX./max(RefImgX)).^x;
        RefImgY=(RefImgY./max(RefImgY)).^x;  
elseif ndims(EM) == 3 % EM - 3D
        RefL{1}=(EM - circshift(EM,[1 0 0]))/BetaVals(1);  % cyclic rotation
        RefL{2}=(EM - circshift(EM,[0 1 0]))/BetaVals(2);  % cyclic rotation
        RefL{3}=(EM - circshift(EM,[0 0 1]))/BetaVals(3);  % cyclic rotation
        RefR{1}=(circshift(EM,[-1 0 0]) - EM)/BetaVals(1);
        RefR{2}=(circshift(EM,[0 -1 0]) - EM)/BetaVals(2);
        RefR{3}=(circshift(EM,[0 0 -1]) - EM)/BetaVals(3);
        RefImgX=sqrt(abssqr(RefL{1})+abssqr(RefL{2})+abssqr(RefL{3}));
        RefImgY=sqrt(abssqr(RefR{1})+abssqr(RefR{2})+abssqr(RefR{3}));
        RefImgX=(RefImgX./max(RefImgX)).^x;
        RefImgY=(RefImgY./max(RefImgY)).^x;               
else % EM - 1D
        RefImgX=sqrt(abssqr((EM - circshift(EM,1))/BetaVals(1)));  
        RefImgY=sqrt(abssqr((circshift(EM,-1) - EM)/BetaVals(1))); 
        RefImgX=(RefImgX./max(RefImgX)).^x;
        RefImgY=(RefImgY./max(RefImgY)).^x;        
end         