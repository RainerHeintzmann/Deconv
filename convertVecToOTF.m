function res=convertVecToOTF(myinput)
global otfrep;  % attention. Will be overwritten
global myim;
global OTFmask;
DataSize=size(myim{1});
% DataLength=prod(DataSize);
PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose

%if length(otfrep) > 1
%    error('Blind PSF deconvolution is currently only implemented for a single PSF');
%end
if isempty(OTFmask) || numel(OTFmask)==0  % This means the vector refers to intensity data
    error('When using blind OTF deconvolution, the global OTFmask needs to be defined');
end

PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose
TotalReadData=0;
for v= 1:numel(otfrep)  % last pattern will be generated from sum-requirement
        if v == 1
            myinput=complex(myinput(1:end/2), myinput(end/2+1:end)); % pack two reals to one complex
        end
        DataToRead=sum(OTFmask{v});
        otfrep{v}=AmpMaskToOTF(OTFmask{v},myinput(1+TotalReadData:TotalReadData+DataToRead));
        TotalReadData=TotalReadData+DataToRead;
        %if ~equalsizes(DataSize,size(otfrep{v}))
        %    error('Datasize unequal to tofill size. This should not happen');
        %end
end

res=otfrep;
