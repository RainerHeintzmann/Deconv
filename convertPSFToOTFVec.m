function out=convertPSFToOTFVec(grad)  % The real-space PSF has to be transformed and converted into a matlab vector
global OTFmask;
global noFFT;

if isempty(noFFT)
    noFFT=0;
end

sumpixels=0;
for v= 1:numel(OTFmask)  % Just to find out the total number of pixels needed to store the data
        sumpixels=sumpixels+sum(OTFmask{v});
end

if exist('zeros_cuda') > 0
    out=zeros_cuda(sumpixels,1,'scomplex'); % allocate the output vector
else
    out=complex(zeros(sumpixels,1),0); % allocate the output vector
end

WrittenData=0;
%currentIlluMaskIdx=1;
currentSumIdx=1;
for v= 1:size(grad,4)  % This loop does the packing
        subgrad=squeeze(grad(:,:,:,v-currentSumIdx));
        csize=size(OTFmask{v});  % Should be size of the rft'ed mask
        %if isa(myillu_mask{v},'cuda')
        csize(2) = (csize(2)-1)*2;  % Calculate back what the corresponding real data size would be
        %end
        if ~equalsizes(size(subgrad),csize) % This means that a border region is used
            subgrad=extract(subgrad,size(OTFmask{v})); % ignore the border region
        end
        
        if noFFT
            transformed=subgrad;
        else
            transformed=rft(subgrad);
            transformed=FixGradRFT(transformed,OTFmask{v});  % Necessary to make the gradient correct for RFTs            
        end
        clear subgrad;        
        toWrite=double(transformed(OTFmask{v}));  % selects only the pixels inside the mask
        WriteSize=numel(toWrite);
        out(1+WrittenData:WrittenData+WriteSize)=toWrite;
        WrittenData=WrittenData+WriteSize;
        clear transformed;
end
out=[real(out);imag(out)]; % unpack complex to two reals
