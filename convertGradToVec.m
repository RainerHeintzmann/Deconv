function out=convertGradToVec(grad)
global ToEstimate;
global myillu_mask;
global myillu_sumcond;

if isempty(ToEstimate) || ToEstimate==0 || isempty(myillu_mask) || numel(myillu_mask) < 1
    out=(reshape(double(grad),[prod(size(grad)) 1]));   % this applies to object as well as unpacked (4D) illumination distributions
    if ~isreal(grad)
        out=[real(out);imag(out)]; % unpack complex to two reals
    end
elseif ToEstimate==1 &&  ~isempty(myillu_mask) && numel(myillu_mask) >= 1
    sumpixels=0;
    currentSumIdx=1;
    for v= 1:numel(myillu_mask)  % Just to find out the total number of pixels needed to store the data
        if v ~= myillu_sumcond{currentSumIdx}
            sumpixels=sumpixels+sum(myillu_mask{v});
        else
        currentSumIdx=currentSumIdx+1;
        end
    end
    %if isa(grad,'cuda')
    if exist('zeros_cuda') > 0
        out=zeros_cuda(sumpixels,1,'scomplex'); % allocate the output vector
    else
        out=complex(zeros(sumpixels,1),0); % allocate the output vector
    end
    %else
        % out=complex(zeros(sumpixels,1,'single'));
    %end
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
                subgrad=extract(subgrad,size(myillu_mask{v})); % ignore the border region
            end
            %if isa(grad,'cuda')
                transformed=rft(subgrad);
            %else
            %    transformed=ft(subgrad);
            %end
            toWrite=double(transformed(myillu_mask{v}));
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
else
    error('Other estimation methos not implemented yet.');
end