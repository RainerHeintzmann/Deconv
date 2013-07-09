function convertVecToIllu(myinput)
global myillu;
global myim;
global myillu_sumcond;
global myillu_mask;
DataSize=size(myim{1});
DataLength=prod(DataSize);
if (isempty(myillu_mask) || numel(myillu_mask)==0) && (numel(myinput) ~= DataLength*(numel(myim)-length(myillu_sumcond)))
    error('The illumination vector is of wrong length. It should not include the last illumination pattern!');
end
asum=0;
sumviews=0;
currentSumCondIdx=1;
PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose
TotalReadData=0;
for v= 1:numel(myim)  % last pattern will be generated from sum-requirement
    if v ~= myillu_sumcond{currentSumCondIdx}
        if isempty(myillu_mask) || numel(myillu_mask)==0
            myillu{v}=dip_image(reshape(myinput(1+DataLength*(v-currentSumCondIdx):DataLength*(v-currentSumCondIdx+1)),PSize),'single');  % The reconstruction can be up to 3D, whereas the data might be 4D
        else  % in this case this is interpreted as complex numers only in the non-zero region of the Fourier-space
            if v == 1
                myinput=complex(myinput(1:end/2), myinput(end/2+1:end)); % pack two reals to one complex
            end
            tmpim=newim(myillu_mask{v},'scomplex');
            DataToRead=sum(myillu_mask{v});  
            tofill=myinput(1+TotalReadData:TotalReadData+DataToRead);   % writes the complex vector into the Fourier-image
            tmpim(myillu_mask{v})=tofill;   % writes the complex vector into the Fourier-image
            TotalReadData=TotalReadData+DataToRead;
            %if isa(tmpim,'cuda')
                tmpim=rift(tmpim);         % go to real space (cuda version)
            %else
            %    tmpim=real(ift(tmpim));    % DipImage version (no cuda). Problem may be that there are double as many variables
            %end
            if ~equalsizes(DataSize,size(tmpim))
                tmpim=extract(tmpim,PSize(1:length(size(tmpim))),[],mean(tmpim));
                error('This should not happen');
            end
            myillu{v}=tmpim;
            clear tmpim;
        end
        asum=asum+myillu{v};
        sumviews=sumviews+1;
    else
        myillu{v}=(sumviews+1)-asum;   % The sum is forced to be constrained
        asum=0;sumviews=0;
        currentSumCondIdx=currentSumCondIdx+1;
    end
end
clear asum;
