function out=convertIlluToVecFourierMask(grad)
global myillu_mask;
global myillu_sumcond;

sumpixels=0;
currentSumIdx=1;
for v= 1:numel(myillu_mask)  % Just to find out the total number of pixels needed to store the data
    if v ~= myillu_sumcond{currentSumIdx}
        sumpixels=sumpixels+sum(myillu_mask{v});
    else
        currentSumIdx=currentSumIdx+1;
    end
end
if exist('zeros_cuda') > 0
    out=zeros_cuda(sumpixels,1,'scomplex'); % allocate the output vector
else
    out=complex(zeros(sumpixels,1),0); % allocate the output vector
end

WrittenData=0;
%currentIlluMaskIdx=1;
currentSumIdx=1;
for v= 1:size(grad,4)+numel(myillu_sumcond)  % This loop does the packing
    if v ~= myillu_sumcond{currentSumIdx};
        subgrad=squeeze(grad(:,:,:,v-currentSumIdx));
        csize=size(myillu_mask{v});  % Should be size of the rft'ed mask
        %if isa(myillu_mask{v},'cuda')
        csize(2) = (csize(2)-1)*2;  % Calculate back what the corresponding real data size would be
        %end
        if ~equalsizes(size(subgrad),csize) % This means that a border region is used
            if length(csize) < 3
                csize(3)=1;
            end
            subgrad=extract(subgrad,csize); % ignore the border region
        end
        transformed=rft(subgrad)*2;  % Aurelie & Rainer to make the gradient correct.  Why ??
        if ndims(transformed)<3
            transformed(:,end)=transformed(:,end)/2; % Aurelie & Rainer to make th egradient correct
        else
            transformed(:,end,:)=transformed(:,end,:)/2; % Aurelie 06.02.2014 for thick slice
        end
        toWrite=double(transformed(myillu_mask{v}));  % selects only the pixels inside the mask
        WriteSize=numel(toWrite);
        out(1+WrittenData:WrittenData+WriteSize)=toWrite;
        WrittenData=WrittenData+WriteSize;
        clear transformed;
        %currentIlluMaskIdx=currentIlluMaskIdx+1;
    else
        currentSumIdx = currentSumIdx+1;
    end
end
out=[real(out);imag(out)]; % unpack complex to two reals
