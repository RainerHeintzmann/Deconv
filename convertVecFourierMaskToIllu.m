% thegrad=convertVecFourierMaskToIllu(myinput) : Converts a Fourier-space kompressed vector of complex values inside the mask into a series of illumination images
% Also images are generated which are implicitely defined via the illumination sum conditions
% thegrad is an empty image of the correct datatype

function res=convertVecFourierMaskToIllu(myinput) 
global myillu;  % will be overwritten
global myim;
global myillu_sumcond;
global myillu_mask;

DataSize=size(myim{1});
if isempty(myillu_mask) || numel(myillu_mask)==0  % This means the vector refers to intensity data
    error('When using Fourier-mask data for illumination, the mask needs to be defined');
end
asum=0;
sumviews=0;
currentSumCondIdx=1;
PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose
TotalReadData=0;
for v= 1:numel(myim)  % last pattern will be generated from sum-requirement
    if v ~= myillu_sumcond{currentSumCondIdx}
        if v == 1
            myinput=complex(myinput(1:end/2), myinput(end/2+1:end)); % pack two reals to one complex
        end
        DataToRead=sum(myillu_mask{v});
        myillu{v}=AmpMaskToIllu(myillu_mask{v},myinput(1+TotalReadData:TotalReadData+DataToRead));
        TotalReadData=TotalReadData+DataToRead;
        if ~equalsizes(DataSize,size(myillu{v}))
            error('Datasize unequal to tofill size. This should not happen');
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
res=myillu;
