% function GenMaskFromIllu(myIllu,MaskThresh,AbberrationTolerance)
% MaskThresh : relative threshold, compared to the maximum pixel intensity  (default: 0.1)
% AbberrationTolerance: a vector with the size in Fourier space to dilate (default [5 5 3])
% writes the generated mask into the global variable myillu_mask;
%
function GenMaskFromIllu(myIllu,MaskThresh,AbberrationTolerance)
if nargin<2
    MaskThresh=0.1;
end
if nargin<2
    AbberrationTolerance=[5 5 3];
end
global myillu_mask;   % confines the variable only to a subspace of illu
myillu_mask=[]; 
% myillu_mask=cell(1,numel(myIllu));
myillu_mask=cell(size(myIllu,1),size(myIllu,2)); % Mask should also be 2*9 Aurelie 27.05.2014
if iscell(myIllu)
%     for v=1:numel(myIllu)
    for si= 1:size(myIllu,1)  % Aurelie 27.05.2014
        for v=1:size(myIllu,2) % Aurelie 27.05.2014
            tomask=abs(ft(myIllu{si,v})); % Aurelie 27.05.2014
%             tomask=abs(ft(myIllu{v}));
            tomask=tomask/max(tomask); % to make it relative
            myillu_mask{si,v}=fft2rft(ifftshift(dilation(tomask > MaskThresh,AbberrationTolerance))); % Rainer's version 21.03.14
        end
    end
else
    tomask=abs(ft(myIllu));
    tomask=tomask/max(tomask); % to make it relative
    myillu_mask=fft2rft(ifftshift(dilation(tomask > MaskThresh,AbberrationTolerance))); % Rainer's version 21.03.14
end