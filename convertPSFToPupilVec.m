function out=convertPSFToPupilVec(grad)  % The real-space PSF has to be transformed and converted into a matlab vector
global OTFmask;
%global noFFT;

%if isempty(noFFT)
%    noFFT=0;
%end

%sumpixels=0;
%for v= 1:numel(OTFmask)  % Just to find out the total number of pixels needed to store the data
%        sumpixels=sumpixels+sum(OTFmask{v});
%end

global PupilInterpolators;  % contains the interpolation coefficients. If this does not exist, it was generated in ParseRegularisation
%global savedATF;
global savedASF;

sumpixels=size(PupilInterpolators.indexList2D,2)*numel(OTFmask);  % A bit of a hack for now

if exist('zeros_cuda') > 0
    out=zeros_cuda(sumpixels,1,'scomplex'); % allocate the output vector
else
    out=complex(zeros(sumpixels,1),0); % allocate the output vector
end


WrittenData=0;
%currentIlluMaskIdx=1;
currentSumIdx=1;
for v= 1:size(grad,4)  % This loop does the packing
        mg=squeeze(grad(:,:,:,v-currentSumIdx));
        if (1)
            if ndims(mg) > 2
                midz=floor(size(mg,3)/2);
                mg=mg(:,:,midz).*savedASF(:,:,midz)*sqrt(size(mg,3));
            else
                mg=mg*savedASF;
            end
            subgrad=conj(ift(conj(mg)));
        else
            mg=mg.*savedASF;  % no factor of 2!
            subgrad=conj(ift(conj(mg)));  % goes to OTF space. The conj is needed due to the Wirtinger Derivative of the Fourier transform
            subgrad=ProjSphere2D(subgrad,PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList);  % estimates the 2D pupil from a given 3D distribution
        end
        % subgrad=2*subgrad.*savedATF;
        toWrite=double(subgrad(PupilInterpolators.indexList2D));  % selects only the pixels inside the mask
        WriteSize=numel(toWrite);
        out(1+WrittenData:WrittenData+WriteSize)=toWrite;
        WrittenData=WrittenData+WriteSize;
        clear subgrad;
end
out=[real(out);imag(out)]; % unpack complex to two reals
