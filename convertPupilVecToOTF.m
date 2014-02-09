function res=convertPupilVecToOTF(myinput)  % converts a vector of the Z-projected pupil function to a 3D OTF
global otfrep;  % attention. Will be overwritten
global myim;
%global OTFmask;  % is supposed to be a 2D pupil mask

%if length(otfrep) > 1
%    error('Blind PSF deconvolution is currently only implemented for a single PSF');
%end
%if isempty(OTFmask) || numel(OTFmask)==0  % This means the vector refers to intensity data
%    error('When using blind Pupil deconvolution, the global OTFmask needs to be given as a 2D image defining the pupil');
%end

global PupilInterpolators;  % contains the interpolation coefficients. If this does not exist, it was generated in ParseRegularisation
global savedATF;
global savedASF;

DataSize=size(myim{1});
PSize=DataSize;
TotalReadData=0;
savedATF=newim(PSize(1:2),'scomplex');  % 2D pupil
atf=newim(PSize,'scomplex');
myinput=complex(myinput(1:end/2), myinput(end/2+1:end)); % pack two reals to one complex

for v= 1:numel(otfrep)  % last pattern will be generated from sum-requirement
        DataToRead=size(PupilInterpolators.indexList2D,2); % sum(OTFmask{v});
        savedATF(PupilInterpolators.indexList2D)=myinput(1+TotalReadData:TotalReadData+DataToRead)*size(myim{1},3);
        atf=FillProjSphere(atf,savedATF,PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList); 
        savedASF=ift(atf);

        otfrep{v}=rft(fftshift(abssqr(savedASF)))/prod(size(myim{1}));
        TotalReadData=TotalReadData+DataToRead;
        %if ~equalsizes(DataSize,size(otfrep{v}))
        %    error('Datasize unequal to tofill size. This should not happen');
        %end
end
if TotalReadData ~= prod(size(myinput))
    error('convertPupilVecToOTF: Data provided is more than expected');
end

res=otfrep;
