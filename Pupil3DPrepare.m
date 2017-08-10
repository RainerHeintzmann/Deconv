% Pupil3DPrepare(imgsize) : optimizes the coefficients for interpolating pupils in 3D and stores preperatory data in the global PupilInterpolators structure
% All the other requires PSF parameters are taken from the PupilInperpolators structure
% PupilInterpolators.Newlambda
% PupilInterpolators.Newpixelsize
% PupilInterpolators.NewRI
% PupilInterpolators.NewNA
% this function also checks whether recomputation is necessary or not.
%
function Pupil3DPrepare(imgsize) % Needs to know the image size. So only possible after accounting for resampling and borders
global PupilInterpolators;
if isempty(PupilInterpolators) || ~isempty(PupilInterpolators.Newlambda) && (~isfield(PupilInterpolators,'psfsize') || norm(imgsize- PupilInterpolators.psfsize) ~= 0 || PupilInterpolators.lambda ~= PupilInterpolators.Newlambda || PupilInterpolators.NA ~= PupilInterpolators.NewNA || PupilInterpolators.RI ~= PupilInterpolators.NewRI || norm(PupilInterpolators.pixelsize - PupilInterpolators.Newpixelsize) ~= 0)
    if exist('disableCuda')
        disableCuda()
    end
    fprintf('Warning: global structure PupilInterpolators does not exist or PSF parameters changed. Recomputing imatrix and pupil factors\n');
    kernelSize=6;  % middle plus 2*4
    % kernelSize=4;  % middle plus 2*2. Seems to be the optimum for 16 slices and 128x128 pixels
    sz=imgsize(3);  % has already been expanded
    Bsize=round(ceil(sz)*0.05);
    imatrix=IterateCoefficients(60,kernelSize,sz,Bsize,500);  % 40 subpixel subdivisions, 2*10+1 kernelsize, 20 pixel bordersize in all directions, 500 iterations
    % imatrix=IterateCoefficients(120,kernelSize,sz,Bsize,500);  % 40 subpixel subdivisions, 2*10+1 kernelsize, 20 pixel bordersize in all directions, 500 iterations
    [indexList2D,fullIndex3D,factorList,aMask]=FillProjSpherePrepare(imgsize,PupilInterpolators.Newlambda,PupilInterpolators.Newpixelsize,PupilInterpolators.NewNA,imatrix,PupilInterpolators.NewRI);
    PupilInterpolators.psfsize=imgsize;
    PupilInterpolators.lambda=PupilInterpolators.Newlambda;
    PupilInterpolators.pixelsize=PupilInterpolators.Newpixelsize;
    PupilInterpolators.NA=PupilInterpolators.NewNA;
    PupilInterpolators.RI=PupilInterpolators.NewRI;
    PupilInterpolators.indexList2D=indexList2D;
    PupilInterpolators.fullIndex3D=fullIndex3D;
    PupilInterpolators.factorList=factorList;
    PupilInterpolators.Mask=aMask;
    ImageParam=struct('Sampling',PupilInterpolators.pixelsize,'Size',imgsize);
    PSFParam=struct('NA',PupilInterpolators.NA,'n',PupilInterpolators.RI,'MinOtf',1.2e-3,'lambdaEm',PupilInterpolators.lambda);
    InputPol=newim([ImageParam.Size(1:2),1,2],'scomplex');
    InputPol(:,:,:,0)=sqrt(0.5);InputPol(:,:,:,1)=i*sqrt(0.5);   % Makes it circular polarisation (or also == random)
    lambdaMedium=PSFParam.lambdaEm/PSFParam.n;
    naAir=(1.0+3*PSFParam.NA/PSFParam.n)/4; % naAir=PSFParam.NA/PSFParam.n;  % This is set to one, as the jinc function will take care of this in a better way
    FPlane=SimLens(InputPol,lambdaMedium,ImageParam.Sampling(1:2),naAir,-1);  % Circular input polarisation, no aplanatic factor, as this is included in the 3D aperture projection method
    PupilInterpolators.Aperture=FPlane .* ft(jincPSF(ImageParam,PSFParam));  % This aperture contains the vectorial effects and the necessary interpolation at the edges of the pupil to avoid stripes
    PupilInterpolators.Aperture=PupilInterpolators.Aperture / abs(PupilInterpolators.Aperture(MidPosX(PupilInterpolators.Aperture),MidPosY(PupilInterpolators.Aperture),:,0));
    if exist('enableCuda')
        enableCuda();
    end
end