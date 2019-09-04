% [RegMat1,RegMat2]=ParseRegularisation(mycells) : Parses the parameter strings into one or two matrices
function [RegMat1,RegMat2,RegMat3]=ParseRegularisation(mycells,toReg)
global DeconvMethod;
global DeconvMask;
global aResampling;
global subSampling;
global myillu;   % here the conversion to cuda could be performed. "useCudaGlobal". Aurelie can you do this?
global myotfs;
global aRecon;
global RefObject;
global myillu_mask;
global myillu_sumcond;
global PupilInterpolators;
global myim;
global realSpaceMultiplier;
global ReadVariance;
global RefImgX;
global RefImgY;
global BwdOTF;

%global ComplexObj;
%ComplexObj=0;
%global IntensityData;
%IntensityData=0;
%global ForcePos;
if nargin<2 
    toReg=0; % meaning object
end

NumMaxPar=35;
RegMat1 = zeros(NumMaxPar,3); % Object updates
RegMat2 = zeros(NumMaxPar,3); % Illumination updates
RegMat3 = zeros(NumMaxPar,3); % PSF updates

DefaultEpsR = 1.0;   % 0.1 is still too small and generates "hot" pixels

switch DeconvMethod
    case 'LeastSqr'
        RegMat1(10,1)=0;
        RegMat2(10,1)=0;
        RegMat3(10,1)=0;
    case 'Poisson'
        RegMat1(10,1)=1;
        RegMat2(10,1)=1;
        RegMat3(10,1)=1;
    case 'WeightedLeastSqr'
        RegMat1(10,1)=2;
        RegMat2(10,1)=2;
        RegMat3(10,1)=2;
    case 'GaussianWithReadnoise'
        RegMat1(10,1)=3;
        RegMat2(10,1)=3;
        RegMat3(10,1)=3;
    case 'Empty'
        RegMat1(10,1)=4;
        RegMat2(10,1)=4;
        RegMat3(10,1)=4;
    otherwise
        fprintf('Unknown Update method: %s\n',DeconvMethod);
        error('For update method only LeastSqr, Poisson and WeigthedLeastSqr are allowed.');
end

if isempty(mycells)
    return
end
if ~iscell(mycells)
    error('For regularisation a cell array has to be provided with the different regularisation names and arguments.')
end
if iscell(mycells{1})  % This means the user wants to parse two or three such arrays
    if numel(mycells) > 3
        error('Use one, two or three cells for regularisation parameters');
    end
    RegMat1=ParseRegularisation(mycells{1},0);
    if numel(mycells)>1 && ~isempty(mycells{2})
        RegMat2=ParseRegularisation(mycells{2},1);
    end
    if numel(mycells)>2 && ~isempty(mycells{3})
        RegMat3=ParseRegularisation(mycells{3},2);
    end
    return
end

for n=1:size(mycells,1)
    switch (mycells{n,1})  % This should be the token
        case 'GS'  % Args are : Lambda
            RegMat1(1,1)=mycells{n,2};
        case 'ER'  % Args are : Lambda.  Used to be called 'AR'
            if numel(mycells{n,2}) ~= 2
                error('Using ER regularisation, please provide two values in the form [lambda eps]. An eps of one is a choice. Eps > 0');
            end
            RegMat1(2,1)=mycells{n,2}(1);
            RegMat1(2,2)=mycells{n,2}(2);
        case 'TV'  % Args are : Lambda, EpsC
            if numel(mycells{n,2}) ~= 2
                error('Using TV regularisation, please provide two values in the form [lambda epsC]. An epsC of Zero means standard TV');
            end
            RegMat1(3,1)=mycells{n,2}(1);
            RegMat1(3,2)=mycells{n,2}(2);
        case 'Kevran'  % Args are : Lambda, EpsC
            RegMat1(32,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) ~= 2
                fprintf('Using Kevran regularisation, please provide two values in the form [lambda epsR]. epsR is assumed to be 0.5 for now.');
                RegMat1(32,2)=0.5;
            else
                RegMat1(32,2)=mycells{n,2}(2);
            end
        case 'NegSqr'  % Args are : Lambda
            RegMat1(4,1)=mycells{n,2};
        case 'GR'  % Args are : Lambda
            RegMat1(5,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(5,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(5,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(5,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(5,3)=0;  % default value: no convolution
            end
        case 'CO'  % Conchello Regularisation. Args are : Lambda
            RegMat1(11,1)=mycells{n,2};
        case 'Reuse'
            RegMat1(6,1)=1;  % Reuse what is in aRecon
        case 'RefObject'  % Reference Object for the Deconv results to be compared to
            RefObject=mycells{n,2};
        case 'StartImg'
            if ~(isa(mycells{n,2},'dip_image') || isa(mycells{n,2},'cuda'))
                error('When submitting a starting image for object or illumination, it has to be a dip_image or cuda type');
            end
            RegMat1(6,1)=1;  % Reuse what is written into aRecon below
            if toReg==0
                aRecon=mycells{n,2};
            elseif toReg==1
                myillu=mycells{n,2};
            elseif toReg==2
                myotfs=mycells{n,2};
            end
        case 'Illumination'
            
            if ~iscell(mycells{n,2})  && (~(isa(mycells{n,2},'dip_image') || isa(mycells{n,2},'cuda')))
                error('When submitting a illumination image for object or illumination, it has to be a dip_image or cuda type');
            end
            myillu=mycells{n,2};
        case 'IlluMask'
            if ~iscell(mycells{n,2}) && (~(isa(mycells{n,2},'dip_image') || isa(mycells{n,2},'cuda')))
                error('When submitting a illumination mask image for object or illumination, it has to be a dip_image or cuda type');
            end
            myillu_mask=mycells{n,2};
        case 'DeconvMask'
            if ~iscell(mycells{n,2}) && (~(isa(mycells{n,2},'dip_image') || isa(mycells{n,2},'cuda')))
                error('When submitting a DeconvMask mask image for residuum, it has to be a dip_image or cuda type');
            end
            DeconvMask=mycells{n,2};
        case 'IlluSums'
            myillu_sumcond=mycells{n,2};
        case 'Complex'
            RegMat1(7,1)=1;  
            %ComplexObj=1;
        case 'IntensityData'
            RegMat1(8,1)=1;              
            %IntensityData=1;
        case 'FTData'  % Data lives in Fourier-space
            RegMat1(21,1)=1;  
        case 'ForcePos'
            RegMat1(9,1)=1;  
            %ForcePos=1;
        case 'Resample'  % reconstructed object is in a different sampling than data
            aResampling=mycells{n,2};
        case 'SubSampling'  % This means the data is subsampled by a fixed factor including aliasing
            subSampling=mycells{n,2};
            if norm(floor(subSampling)-subSampling) ~=0
                error('Only Integer values are allowed for "subSampling"');
            end
        case 'Bg'  % include an offset intensity into the model
            RegMat1(12,1)=mycells{n,2};
        case 'Show' 
            RegMat1(13,1)=1;
        case 'ProjPupil'   % This means the only the projected pupil is iterated rather that the full 3d amplitude. 
            RegMat1(14,1)=1;
            if (size(mycells,2) < 3) || (numel(mycells{n,2}) ~= 3) || (numel(mycells{n,3}) ~= 3)
                error('Error using the argument ''ProjPupil''. You need to state wavelength and NA as follows: {''ProjPupil'', [lambda, NA, n],[pixelsizeX pixelsizeY pixelsizeZ]}');
            end
            PupilInterpolators.Newlambda=mycells{n,2}(1);  % Wavelength
            PupilInterpolators.NewNA=mycells{n,2}(2);  % NA
            PupilInterpolators.NewRI=mycells{n,2}(3);  % refractive index
            PupilInterpolators.Newpixelsize=mycells{n,3};
            PupilInterpolators.FullVectorial=0;  % Scalar approach
            if mycells{n,2}(1)== 2  % use full vectorial theory
                PupilInterpolators.FullVectorial=1;
            end
        case 'ForcePhase'
            RegMat1(15,1)=1;  
            %ForcePos=1;
        case 'NormMeasSum'
            RegMat1(16,1)=1;
        case 'NormMeasSumSqr'
            RegMat1(17,1)=1;
        case 'NormFac'  
            RegMat1(18,1)=mycells{n,2}(1);
        case 'MaxTestDim'  % Number of dimensions to test for gradient tests
            RegMat1(19,1)=mycells{n,2}(1);
        case 'RealSpaceMask'  % used for light sheet microscopy in the PSF regularisation field. If only a number is supplied, this is interpreted as the width of the Gaussian
            if prod(size(mycells{n,2}(1))) == 1
                RegMat1(20,1)=toReg+1;  % just serves as an index
                realSpaceMultiplier{toReg+1}=exp(-zz(size(myim{1})).^2./(log(0.5)*mycells{n,2}(1)).^2);
            else
                RegMat1(20,1)=toReg+1;
                realSpaceMultiplier{toReg+1}=mycells{n,2}(1);    
            end
        case 'LAP'  % Args are : Lambda
            RegMat1(22,1)=mycells{n,2}(1);
        case 'GRLapGrad6'  % Args are : Lambda
            RegMat1(23,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(23,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(23,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(23,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(23,3)=0;  % default value: no convolution
            end 
            
        case 'GRLapDiv6'  % Args are : Lambda
            RegMat1(29,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(29,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(29,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(29,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(29,3)=0;  % default value: no convolution
            end
        case 'GRLapDiv6stabil'  % Args are : LambdaRegularize
            RegMat1(30,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(30,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(30,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(30,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(30,3)=0;  % default value: no convolution
            end
        case 'GRLapDivReg2'  % Args are : LambdaRegularize
            RegMat1(31,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(31,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(31,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(31,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(31,3)=0;  % default value: no convolution
            end
            
        case 'GRLapGradReg'  % Args are : Lambda
            RegMat1(24,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(24,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(24,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(24,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(24,3)=0;  % default value: no convolution
            end  
        case 'GRLapDivReg'  % Args are : Lambda
            RegMat1(28,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(28,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(28,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(28,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(28,3)=0;  % default value: no convolution
            end
        case 'GRCentral'  % Args are : Lambda
            RegMat1(25,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(25,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(25,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(25,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(25,3)=0;  % default value: no convolution
            end  
       case 'GRLap6'  % Args are : Lambda
            RegMat1(26,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(26,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(26,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(26,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(26,3)=0;  % default value: no convolution
            end
        case 'Lap27'  % Args are : Lambda
            RegMat1(27,1)=mycells{n,2}(1);
            if numel(mycells{n,2}) > 1
                RegMat1(27,2)=mycells{n,2}(2);  % optional modified Good's roughness
            else
                RegMat1(27,2)=DefaultEpsR; % This number is too low and causes hot pixels: 1e-4;  % default value for nominator regularisation
            end
            if numel(mycells{n,2}) > 2
                RegMat1(27,3)=mycells{n,2}(3);  % optional flag for convolution in the nominator
            else
                RegMat1(27,3)=0;  % default value: no convolution
            end 
        case 'CLE_GS'  % Args are : Lambda
            RegMat1(33,1)=mycells{n,2}{1};
            RefImgX=mycells{n,2}{2};  % Reference Img X
            RefImgY=mycells{n,2}{3};  % Reference Img Y
        case 'ReadVariance'
            ReadVariance=mycells{n,2}(1);
        case 'NoPSF'
            RegMat1(34,1)=1;  % will disable both convolution operations in the forward and backward models
        case 'WienerBwd'
            RegMat1(35,1)=mycells{n,2};  % will replace the backward pass with a Wiener-filtered version (Shroff paper 2019).
        otherwise
            fprintf('Unknown Flag: %s\n',mycells{n,1});
            error('For regularisation only TV, ER, GR, CO, Kevran, Complex, IntensityData, ForcePos, ForcePhase, NormMeasSum, NormMeasSumSqr, NormFac, MaxTestDim, NegSqr, Reuse, Resample, Bg, ProjPupil, Illumination, IlluMask, FTData and StartImg are allowed');
    end
end
