function out=convertIlluToVecFourierMask(grad)
global myillu_mask;
global myillu_sumcond;

sumpixels=0;
currentSumIdx=1;
for si= 1:size(myillu_mask,1)  % sub illuminations
    for v= 1:size(myillu_mask,2)  % Just to find out the total number of pixels needed to store the data
%         if v ~= myillu_sumcond{currentSumIdx}
        if v ~= myillu_sumcond{si,currentSumIdx} % Aurelie 23.06.2014. Theoretically, the number of sub illuminations should correspond to the first dimension of the illumination sum condition as well (this probably requires a sanity check)
            sumpixels=sumpixels+sum(myillu_mask{si,v});
        else
            currentSumIdx=currentSumIdx+1;
        end
    end
    currentSumIdx=1; % Aurelie 23.06.2014. For the loop on si, start from 1 again
end
if exist('zeros_cuda') > 0
    out=zeros_cuda(sumpixels,1,'scomplex'); % allocate the output vector
else
    out=complex(zeros(sumpixels,1),0); % allocate the output vector
end

WrittenData=0;
%currentIlluMaskIdx=1;
currentSumIdx=1;
if ndims(grad) < 4
    gradsize=size(grad,3);
else
    gradsize=size(grad,4);
end
% for v= 1:gradsize+numel(myillu_sumcond)  % This loop does the packing
for v= 1:gradsize+size(myillu_sumcond,2)  % This loop does the packing. Aurelie 23.06.2014
%     if v ~= myillu_sumcond{currentSumIdx};
    if v ~= myillu_sumcond{1,currentSumIdx}; % Aurelie 23.06.2014. Not sure whether we should do a loop on si=size(myillu_sumcond,1)
        if ndims(grad) < 4 %Aurelie 10.03.2014. 2D blind-SIM with mask
            subgrad=squeeze(grad(:,:,v-currentSumIdx));
        else
            subgrad=squeeze(grad(:,:,:,v-currentSumIdx));
        end
%         csize=size(myillu_mask{si,v});  % Should be size of the rft'ed mask
        csize=size(myillu_mask{1,1});  % Should be size of the rft'ed mask. Aurelie 27.05.2014: this is just to get the size, and v may be bigger than the number of elements in the mask
        %if isa(myillu_mask{v},'cuda')
        csize(2) = (csize(2)-1)*2;  % Calculate back what the corresponding real data size would be
        %end
        if ndims(csize) < 3 && ndims(subgrad) > 2  % For the case of a two-D mask, which is interpreted as the central slice in Fourier-space
            csize(3) = size(subgrad,3);
        end
        if ~equalsizes(size(subgrad),csize) % This means that a border region is used
            if ndims(subgrad) > 2 && length(csize) < 3
                csize(3)=1;
            end
            subgrad=extract(subgrad,csize); % ignore the border region
        end
        transformed=rft(subgrad);  
        transformed=FixGradRFT(transformed,myillu_mask{1,v});  % Necessary to make the gradient correct for RFTs            
        for si= 1:size(myillu_mask,1)  % sub illuminations
            toWrite=double(transformed(myillu_mask{si,v}));  % selects only the pixels inside the mask. Remark from Aurelie: transformed and myillu_mask have to be of the same type (both cuda or both not)
            WriteSize=numel(toWrite);
            out(1+WrittenData:WrittenData+WriteSize)=toWrite;
            WrittenData=WrittenData+WriteSize;
        end
        clear transformed;
        %currentIlluMaskIdx=currentIlluMaskIdx+1;
    else
        currentSumIdx = currentSumIdx+1;
    end
end
out=[real(out);imag(out)]; % unpack complex to two reals
