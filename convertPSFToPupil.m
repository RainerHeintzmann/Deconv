function out=convertPSFToPupil(grad,ApplyMask)  % The real-space PSF has to be transformed and converted into a matlab vector

global PupilInterpolators;  % contains the interpolation coefficients. If this does not exist, it was generated in Pupil3DPrepare
global savedASF;
% global savedATF;
global realSpaceMultiplier;

currentSumIdx=1;
grad=expanddim(grad,4); % just in case this is only a single view
out={};
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
            out{n}=squeeze(sum(subgrad.*conj(PupilInterpolators.Aperture),[],4));
        end
end
out = cat(4,out);
