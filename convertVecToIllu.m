% function res=convertVecToIllu(myinput) : Converts a vector of real-space illumination data 

function res=convertVecToIllu(myinput)

global myillu; % is overwritten!! To save memory
global myim;
global myillu_sumcond;
global myillu_mask;

DataSize=size(myillu{1});
DataLength=prod(DataSize);
if ~(isempty(myillu_mask) || numel(myillu_mask)==0 ) % This means the vector refers to intensity data
    error('A Fourier-mask is defined even though the vector represents real-space data. Please clear the global variable myillu_mask.');
end

if (numel(myinput) ~= DataLength*(numel(myim)-length(myillu_sumcond)))
    error('The illumination vector is of wrong length. It should not include the illumination patterns as given by the myillu_sumcond!');
end
sumviews=0;
currentSumCondIdx=1;
PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose
for v= 1:numel(myim)  % last pattern will be generated from sum-requirement
    if isempty(myillu_sumcond) || v ~= myillu_sumcond{currentSumCondIdx}
        myillu{v}=dip_image(reshape(myinput(1+DataLength*(v-currentSumCondIdx):DataLength*(v-currentSumCondIdx+1)),PSize),'single');  % The reconstruction can be up to 3D, whereas the data might be 4D
        sumviews=sumviews+1;
    else
        sumviews=0;  % the sum condition is applied in a later stage, since there may be other things happening before (e.g. fft or ForcePos).
        currentSumCondIdx=currentSumCondIdx+1;
    end
end

res=myillu;