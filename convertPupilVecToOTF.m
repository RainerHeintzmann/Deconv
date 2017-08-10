% res=convertPupilVecToOTF(myinput,ApplyMask)  % converts a vector of the Z-projected pupil function to a 3D OTF
% otfrep is assigned and caculated from the incoherent PSF to the global variable and returned by this function
% savedATF and savedASF  are global variables that are assigned to store the information for later calculation of the gradient
% in the function convertPSFToPupilVec

function res=convertPupilVecToOTF(myinput,ApplyMask)  % converts a vector of the Z-projected pupil function to a 3D OTF
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
global realSpaceMultiplier;

DataSize=size(myim{1});
PSize=DataSize;
TotalReadData=0;
savedATF=newim(PSize(1:2),'scomplex');  % 2D pupil
myinput=complex(myinput(1:end/2), myinput(end/2+1:end)); % pack two reals to one complex

for v= 1:numel(otfrep)  % last pattern will be generated from sum-requirement
        DataToRead=size(PupilInterpolators.indexList2D,2); % sum(OTFmask{v});
        %  savedATF(PupilInterpolators.indexList2D)=myinput(1+TotalReadData:TotalReadData+DataToRead)*size(myim{1},3);
        tmp=myinput(1+TotalReadData:TotalReadData+DataToRead)*size(myim{1},3); % other version is not yet Cuda-Compatible
        savedATF(PupilInterpolators.indexList2D)=tmp;
        clear tmp;
        atf=newim(PSize,'scomplex');
        savedATF=savedATF.*PupilInterpolators.Aperture;
        % Now the forward model of the ASF and then PSF is calculated either in the scalar or the vectorial way:
        if ~PupilInterpolators.FullVectorial
            atf=FillProjSphere(atf,savedATF,PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList); 
            savedASF=ift(atf);
            mypsf = abssqr(savedASF);
        else
            atfX=FillProjSphere(atf,squeeze(savedATF(:,:,:,0)),PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList);
            atfY=FillProjSphere(atf,squeeze(savedATF(:,:,:,1)),PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList);
            atfZ=FillProjSphere(atf,squeeze(savedATF(:,:,:,2)),PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList);
            amp3d=cat(4,ift(atfX),ift(atfY),ift(atfZ));  % 4th dimension stores the X, Y and Z field vector
            %amp3d=ift(atfX);  % 4th dimension stores the X, Y and Z field vector
            savedASF=amp3d;
            mypsf=squeeze(sum(abssqr(amp3d),[],4));
            %mypsf=abssqr(amp3d);
        end
        if ApplyMask
            mypsf=mypsf .* realSpaceMultiplier{ApplyMask};
        end
        
        otfrep{v}=rft(fftshift(mypsf))/prod(size(myim{1}));
        clear mypsf;
        TotalReadData=TotalReadData+DataToRead;
        %if ~equalsizes(DataSize,size(otfrep{v}))
        %    error('Datasize unequal to tofill size. This should not happen');
        %end
end
if TotalReadData ~= prod(size(myinput))
    error('convertPupilVecToOTF: Data provided is more than expected');
end

res=otfrep;
