function out=convertPSFToPupilPhaseVec(grad,ApplyMask,isGradient)  % The real-space PSF has to be transformed and converted into a matlab vector
global OTFmask;
global PupilInterpolators;  % contains the interpolation coefficients. If this does not exist, it was generated in Pupil3DPrepare
global savedASF;
global realSpaceMultiplier;
global savedPhaseInput
global AmpFactor;

sumpixels=size(PupilInterpolators.indexList2D,2)*numel(OTFmask);  % A bit of a hack for now

if exist('zeros_cuda','file') > 0
    out=zeros_cuda(sumpixels,1,'single'); % allocate the output vector
else
    out=zeros(sumpixels,1); % allocate the output vector
end

if nargin < 3
    isGradient=1;
end

WrittenData=0;
%currentIlluMaskIdx=1;
currentSumIdx=1;
grad=expanddim(grad,4); % just in case this is only a single view
for v= 1:size(grad,4)  % This loop does the packing
        mg=ifftshift(squeeze(grad(:,:,:,v-currentSumIdx))) / sqrt(prod(size(grad)));  % ifftshift
        if ApplyMask
            mg=mg .* realSpaceMultiplier{ApplyMask};
        end

        if (0)
            if ndims(mg) > 2
                midz=0; % floor(size(mg,3)/2);
                mg=2*mg(:,:,midz).*fftshift(savedASF(:,:,floor(size(mg,3)/2)))*sqrt(size(mg,3));
                %subgrad=ift(mg(:,:,midz)) .*savedATF;
            else
                mg=mg*savedASF;
            end
            %subgrad=conj(ift(conj(mg)))/ sqrt(prod(size(mg)));
            subgrad=ft(fftshift(mg))/ sqrt(prod(size(mg)));
        else
            mg=2*repmat(mg,[1 1 1 size(savedASF,4)]).*savedASF;  % can this be somehow accelerated?
            %subgrad=conj(ift(conj(mg)));  % goes to OTF space. The conj is needed due to the Wirtinger Derivative of the Fourier transform
            subgrad=ft3d(mg) * (size(mg,3)) / sqrt(prod(size3d(mg)));  % goes to OTF space. The conj is needed due to the Wirtinger Derivative of the Fourier transform
            subgrad=ProjSphere2D(subgrad,PupilInterpolators.indexList2D,PupilInterpolators.fullIndex3D,PupilInterpolators.factorList,1);  % estimates the 2D pupil from a given 3D distribution
            subgrad=squeeze(sum(subgrad.*conj(PupilInterpolators.Aperture),[],4));
        end
        subgrad=AmpFactor*double(subgrad(PupilInterpolators.indexList2D))';  % not sure, if the AmpFactor is correct here...
        if isGradient
            % subgrad = -sin(savedPhaseInput{v}).*real(subgrad) + cos(savedPhaseInput{v}).*imag(subgrad);
            subgrad = sin(savedPhaseInput{v}).*real(subgrad) - cos(savedPhaseInput{v}).*imag(subgrad);
        else % inverse
            subgrad = angle(subgrad);
        end
        
        % subgrad=2*subgrad.*savedATF;
        toWrite=subgrad;  % selects only the pixels inside the mask
        WriteSize=numel(toWrite);
        out(1+WrittenData:WrittenData+WriteSize)=toWrite;
        WrittenData=WrittenData+WriteSize;
        clear subgrad;
end
% output is already real
