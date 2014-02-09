function out=convertIlluToVec(grad)
global myillu_sumcond;

s=size(grad);
sumpixels=prod(s);  % Will not contain the implicit illumination data

if isreal(grad)
    if exist('zeros_cuda','file') > 0
        out=zeros_cuda(sumpixels,1,'single'); % allocate the output vector
    else
        out=zeros(sumpixels,1); % allocate the output vector
    end
else
    if exist('zeros_cuda','file') > 0
        out=zeros_cuda(sumpixels,1,'scomplex'); % allocate the output vector
    else
        out=complex(zeros(sumpixels,1),0); % allocate the output vector
    end
end
WrittenData=0;
%currentIlluMaskIdx=1;
currentSumIdx=1;
for v= 1:size(grad,4)+numel(myillu_sumcond)  % This loop does the packing
    if v ~= myillu_sumcond{currentSumIdx};
        subgrad=squeeze(grad(:,:,:,v-currentSumIdx));
        toWrite=double(subgrad);  % selects only the pixels inside the mask
        WriteSize=numel(toWrite);
        out(1+WrittenData:WrittenData+WriteSize)=toWrite;
        WrittenData=WrittenData+WriteSize;
    else
        currentSumIdx = currentSumIdx+1;
    end
end
if ~isreal(grad)
    out=[real(out);imag(out)]; % unpack complex to two reals
end