% function res=convertVecToIllu(myinput) : Converts a vector of real-space illumination data 

function res=convertVecToIllu(myinput) 

global myillu; % is overwritten!! To save memory
global myim;
global myillu_sumcond;
global myillu_mask;
global savedInput;
global RegularisationParameters;

DataSize=size(myillu{1});
DataLength=prod(DataSize);
if ~(isempty(myillu_mask) || numel(myillu_mask)==0 ) % This means the vector refers to intensity data
    error('A Fourier-mask is defiend even though the vector represents real-space data. Please clear the global variable myillu_mask.');
end

if (numel(myinput) ~= DataLength*(numel(myim)-length(myillu_sumcond)))
    error('The illumination vector is of wrong length. It should not include the illumination patterns as given by the myillu_sumcond!');
end
asum=0;
sumviews=0;
currentSumCondIdx=1;
PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose
for v= 1:numel(myim)  % last pattern will be generated from sum-requirement
    if isempty(myillu_sumcond) || v ~= myillu_sumcond{currentSumCondIdx}
        myillu{v}=dip_image(reshape(myinput(1+DataLength*(v-currentSumCondIdx):DataLength*(v-currentSumCondIdx+1)),PSize),'single');  % The reconstruction can be up to 3D, whereas the data might be 4D
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