% AssignFunctions(RegulisierungsParameters,ToEstimate) : Assigns the helper functions FwdModel, BwdModel and CalcResiduum, depending on the model desctibed by ToEstimate
% ToEstimate = 0 : Object only estimation
% ToEstimate = 1 : Object estimation in a model including illumination
% ToEstimate = 2 : Illumination estimation
% ToEstimate = 3 : PSF estimation
function AssignFunctions(RegularisationParameters,ToEstimate)

global FwdModel; % A function (pointer) which is assigned outside. Options are FwdObjConvPSF() FwdObjConvASF() FwdObjIlluConvPSF() and FwdObjIlluConvASF()
global BwdModel;  % This performs the convolution of the residuum with the psf. Options are BwdModel() are BwdResidObjConvPSF() BwdResidObjConvASF() BwdResidObjIlluConvPSF() BwdResidObjIlluConvASF() BwdResidIlluObjConvPSF() BwdResidIlluObjConvASF()
global BwdExtraModel;  % Here changes to the Bwd projection are applied to, such as Sqr or Phase
global CalcResiduum; % A function (pointer) which is assigned outside. Options are ResidPoisson(), ResidLeastSrq(), ResidWeightedLeastSqr()
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF
global ConvertGradToVec; % Contains the routine to convert the model (e.g. the gradient of the object) into the vector to iterate
global ConvertModelToVec; % Contains the routine to convert the model (e.g. the gradient of the object) into the vector to iterate
global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
global myim; % This is the measured data. It is needed to see, if it is of type complex, in which case even real reconstructions have to generate full complex data.
global myFTim; % This is the measured data. It is needed to see, if it is of type complex, in which case even real reconstructions have to generate full complex data.
global useFTComparison; % This is set to one in Ptychography (fast) updates, otherwise zero
global otfrep;  % This is needed to compare the Z-size, to see whether the thick slice (single plane) speedup trick can be used
global ComplexPSF; % will be set, if PSF is complex valued
global myillu_mask; % to decide which Vector conversion routine to use
global my_sumcond;
global myillu_sumcond;

%if RegulisierungsParameters(7,1) % case 'Complex'
%if RegulisierungsParameters(8,1) % case 'IntensityData'
%if RegulisierungsParameters(9,1) % case 'ForcePos'

if RegularisationParameters(10,1) == 0
    CalcResiduum=@ResidLeastSqr;
elseif RegularisationParameters(10,1) == 1
    CalcResiduum=@ResidPoisson;
elseif RegularisationParameters(10,1) == 2
    CalcResiduum=@ResidWeightedLeastSqr;
elseif RegularisationParameters(10,1) == 3
    CalcResiduum=@ResidGaussianWithReadnoise;
elseif RegularisationParameters(10,1) == 4
    CalcResiduum=@ResidZero;
else
    error('Unknown residuum calculation method');
end
% if RegularisationParameters(21,1)  % 'FTData'
%     if ~RegularisationParameters(7,1) % 'complex'
%         error('for FTData, a complex valued object should be reconstructed: Choose flag complex')
%     end
% end

FwdModel={};BwdModel={};ConvertInputToModel={};ConvertModelToVec={};ConvertGradToVec={};
my_sumcond={};
if ToEstimate==0  % Object is estimated, illumination and psf are assumed to be known
    AssignToGlobal=@AssignToObject;
    if ~isreal(myim{1}) || ComplexPSF % data is of type complex or the "IntensityData" flag is activated
        if RegularisationParameters(21,1)  % 'FTData'
            FwdModel=@FwdObjFT;
            BwdModel=@BwdResidObjFT;
        else
            FwdModel=@FwdObjConvASF;
            BwdModel=@BwdResidObjConvASF;
        end
        if RegularisationParameters(8,1)
		if ~isreal(myim{1})
			error('IntensityData flag was used with complex valued measurements.');
		end
           FwdModel = @(aRecon,ftRecon,myIllum,myOtf,norm3D)FwdCalcInt(aRecon,ftRecon,myIllum,myOtf,norm3D, FwdModel);
           BwdModel = @(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)BwdCalcInt(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D, BwdModel);
        end
        if RegularisationParameters(7,1) % case 'Complex'
            ConvertInputToModel=@convertVecToCpxObj;  % This is a packed complex vector
            ConvertModelToVec=@convertCpxObjToVec;
            if RegularisationParameters(15,1)
                error(['Do not use the flag "Complex" for reconstruction of phase-only (ForcePhase-) objects. '...
                    'You have to provide a complex valued asf for phase-only reconstruction. '...
                    'If your data is real valued due to intensity measurement, please set the "IntensityData" flag.']);
            end
        else
            if RegularisationParameters(15,1) % ForcePhase, Forces object to be a phase only object
                ConvertInputToModel=@convertVecToObj;  
                ConvertModelToVec=@convertObjToVec;
            else
                ConvertInputToModel=@convertRVecToCpxObj;  % This is only a real valued vector but expanded to complex.
                ConvertModelToVec=@convertCpxObjToRVec;
            end
        end
    else  % data is real valued. Thus PSF is also real
        if RegularisationParameters(7,1) % 'Complex' reconstruction but real valued data
            if RegularisationParameters(15,1)
                error(['Do not use the flag "Complex" for reconstruction of phase-only (ForcePhase-) objects. '...
                    'You have to provide a complex valued asf for phase-only reconstruction. '...
                    'If your data is real valued due to intensity measurement, please set the "IntensityData" flag.']);
            end
            if RegularisationParameters(8,1) % 'Intensity' data. In this case a complex reconstruction can be attempted
                FwdModel=@FwdObjConvASFSqr;
                BwdModel=@BwdResidObjConvASFSqr;
                ConvertInputToModel=@convertVecToCpxObj;  % This is a packed complex vector
                ConvertModelToVec=@convertCpxObjToVec;
            else
                error('Flag Complex does not make sense for real valued data, if the Intensity flag is not set.');
            end
        else % real reconstruction and real data
            ConvertInputToModel=@convertVecToObj;
            if ndims(otfrep{1}) > 2 && size(myim{1},3) == 1 && size(otfrep{1},3) ~= size(myim{1},3)
                FwdModel=@FwdObjConvPSF_Slice;
                BwdModel=@BwdResidObjConvPSF_Slice;
            else
               FwdModel=@FwdObjConvPSF;
               BwdModel=@BwdResidObjConvPSF;
            end
            ConvertModelToVec=@convertObjToVec;
        end
    end
    if RegularisationParameters(34,1)  % 'NoPSF'
        FwdModel = @FwdIdentity;
        BwdModel = @BwdIdentity;
    end
elseif ToEstimate==1  % Object is estimated, illumination is assumed known, but spatially variing
     AssignToGlobal=@AssignToObject;
    if ~isreal(myim{1}) || ComplexPSF || RegularisationParameters(21,1) % data is of type complex or the "IntensityData" flag is activated
        if RegularisationParameters(21,1)  % 'FTData'
            FwdModel=@FwdObjIlluFT;
            BwdModel=@BwdResidObjIlluFT;
        else
            FwdModel=@FwdObjIlluConvASF;
            BwdModel=@BwdResidObjIlluConvASF;
        end
        if RegularisationParameters(8,1)
		if ~isreal(myim{1})
			error('IntensityData flag was used with complex valued measurements.');
		end
           FwdModel = @(aRecon,ftRecon,myIllum,myOtf,norm3D)FwdCalcInt(aRecon,ftRecon,myIllum,myOtf,norm3D, FwdModel);
           BwdModel = @(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D)BwdCalcInt(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D, BwdModel);
        end
        if RegularisationParameters(7,1) % case 'Complex'
            ConvertInputToModel=@convertVecToCpxObj;  % This is a packed complex vector
            ConvertModelToVec=@convertCpxObjToVec;
            if RegularisationParameters(15,1)
                error(['Do not use the flag "Complex" for reconstruction of phase-only (ForcePhase-) objects. '...
                    'You have to provide a complex valued asf for phase-only reconstruction. '...
                    'If your data is real valued due to intensity measurement, please set the "IntensityData" flag.']);
            end
        else
            if RegularisationParameters(15,1) % ForcePhase, Forces object to be a phase only object
                ConvertInputToModel=@convertVecToObj;  
                ConvertModelToVec=@convertObjToVec;
            else
                ConvertInputToModel=@convertRVecToCpxObj;  % This is only a real valued vector but expanded to complex.
                ConvertModelToVec=@convertCpxObjToRVec;
            end
        end
    else  % data is real valued. Thus PSF is also real
        if RegularisationParameters(7,1) % 'Complex' reconstruction but real valued data
            if RegularisationParameters(15,1)
                error(['Do not use the flag "Complex" for reconstruction of phase-only (ForcePhase-) objects. '...
                    'You have to provide a complex valued asf for phase-only reconstruction. '...
                    'If your data is real valued due to intensity measurement, please set the "IntensityData" flag.']);
            end
            if ~RegularisationParameters(21,1)  % When no PSFs are used (FTData) a Ptychography reconstruction should work
                if RegularisationParameters(8,1) % 'Intensity' data. In this case a complex reconstruction can be attempted
                    error('Intensity Data and real PSF cannot be combined with illumination at the moment.');
                else
                    error('Flag Complex does not make sense for real valued data, if the Intensity flag is not set.');
                end
            end
        else % real reconstruction and real data
            ConvertInputToModel=@convertVecToObj;
            if ndims(otfrep{1}) > 2 && size(myim{1},3) == 1 && size(otfrep{1},3) ~= size(myim{1},3)
                FwdModel=@FwdObjIlluConvPSF_Slice;
                BwdModel=@BwdResidObjIlluConvPSF_Slice;
            else
                FwdModel=@FwdObjIlluConvPSF;
                BwdModel=@BwdResidObjIlluConvPSF;
            end
            ConvertModelToVec=@convertObjToVec;
        end
    end
    %%%%%%
%     AssignToGlobal=@AssignToObject;
%     
%     if ~isreal(myim{1}) || RegularisationParameters(7,1) % case 'Complex'
%         FwdModel=@FwdObjIlluConvASF;
%         BwdModel=@BwdResidObjIlluConvASF;
%     else
%         if ndims(otfrep{1}) > 2 && size(myim{1},3) == 1 && size(otfrep{1},3) ~= size(myim{1},3)  % Can the slice trick be used?
%             FwdModel=@FwdObjIlluConvPSF_Slice;
%             BwdModel=@BwdResidObjIlluConvPSF_Slice;
%         else
%             FwdModel=@FwdObjIlluConvPSF;
%             BwdModel=@BwdResidObjIlluConvPSF;
%         end
%     end
elseif ToEstimate==2  % object and psf are assumed to be known, and illumination is estimated
    my_sumcond=myillu_sumcond;   % to ensure the sumcondition is used
    AssignToGlobal=@AssignToIllu;
    if (isempty(myillu_mask) || numel(myillu_mask)==0 )
        ConvertInputToModel=@convertVecToIllu;
        ConvertModelToVec=@convertIlluToVec;
    else
        ConvertInputToModel=@convertVecFourierMaskToIllu;
        ConvertModelToVec=@convertIlluToVecFourierMask;
    end
    if RegularisationParameters(7,1) % case 'Complex'
        FwdModel=@FwdObjIlluConvASF;
        BwdModel=@BwdResidIlluObjConvASF;
    else
        if ndims(otfrep{1}) > 2 && size(myim{1},3) == 1 && size(otfrep{1},3) ~= size(myim{1},3)
            FwdModel=@FwdObjIlluConvPSF_Slice;
            BwdModel=@BwdResidIlluObjConvPSF_Slice;
        else
            FwdModel=@FwdObjIlluConvPSF;
            BwdModel=@BwdResidIlluObjConvPSF;
        end
    end
elseif ToEstimate==3  % object and illumination are assumed to be known, and otf is estimated
    if RegularisationParameters(14,1) % case 'ProjPupil' where only the 2D pupil is estimated
        AssignToGlobal=@AssignToOTF;
        ConvertInputToModel=@(vec)convertPupilVecToOTF(vec,RegularisationParameters(20,1));  % the Regularistion parameter determines whether to apply a multiplication with a global (Gaussian) mask
        FwdModel=@FwdObjConvPSF;   % identical for object and PSF (already converted to OTF)
        BwdModel=@BwdResidPSFConvObj;  % has to convolve with aRecon
        ConvertModelToVec=@(vec)convertPSFToPupilVec(vec,RegularisationParameters(20,1));  % this conversion is identical for object or PSF
    else
        AssignToGlobal=@AssignToOTF;
        ConvertInputToModel=@convertVecToOTF;  % fills a volume in Fourier-space
        FwdModel=@FwdObjConvPSF;   % identical for object and PSF (already converted to OTF)
        BwdModel=@BwdResidPSFConvObj;  % has to convolve with aRecon
        ConvertModelToVec=@convertPSFToOTFVec;  % this conversion is identical for object or PSF
    end
end

    
if isempty(ConvertGradToVec);
    ConvertGradToVec=ConvertModelToVec;    
end

BwdExtraModel=@(grad,viewNum)NoChange(grad,viewNum);

if RegularisationParameters(15,1) % ForcePhase, Forces object to be a phase only object
    ConvertInputToModel=@(vec)ApplyPhaseModel(vec,ConvertInputToModel);  % fixes the second parameter is the current model
    % ConvertGradToVec=@(grad)ApplyPhaseModelGrad(grad,ConvertGradToVec);
    % BwdModel=@(grad)ApplyPhaseModelBwd(grad,viewNum,BwdModel);
    BwdExtraModel=@(grad,viewNum)ApplyPhaseModelBwd(grad,viewNum,BwdExtraModel);
    ConvertModelToVec=@(model)ApplyInversePhaseModel(model,ConvertModelToVec);
    if RegularisationParameters(9,1)
        error('Please decide whether you want to reconstruct a phase-only (ForcePhase-) or a nonnegativ (ForcePos-) object');
    end
elseif RegularisationParameters(9,1) % ForcePos, Forces object to be positive by estimating only its squareroot
    ConvertInputToModel=@(vec)ApplySqrModel(vec,ConvertInputToModel);  % fixes the second parameter is the current model
    %ConvertGradToVec=@(grad)ApplySqrModelGrad(grad,ConvertGradToVec);
    BwdExtraModel=@(grad,viewNum)ApplySqrModelBwd(grad,viewNum,BwdExtraModel);
    ConvertModelToVec=@(model)ApplyInverseSqrModel(model,ConvertModelToVec);
end
% if RegularisationParameters(20,1)~=0 % RealSpaceMask  : E.g. used for light sheet blind PSF deconvolution to multiply the PSF with a known Z-Gaussian    
%     ConvertInputToModel=@(vec)ApplyMaskModel(vec,ConvertInputToModel,RegularisationParameters(20,1));  
%     BwdExtraModel=@(grad,viewNum)ApplyMaskModelBwd(grad,viewNum,BwdExtraModel,RegularisationParameters(20,1));
%     ConvertModelToVec=@(model)ApplyInverseMaskModel(model,ConvertModelToVec,RegularisationParameters(20,1));
% end

if RegularisationParameters(16,1)  % NormMeasSum, normalizes the result of the FwdModel befor calculating the error and residuum    
    if ~isreal(myim{1})
        error('NormMeasSum can only be applied to real valued measurements. Choose NormMeasSumSqr for complex valued data');
    end
    FwdModel=@(aRecon2,ftRecons,myIllum,myOtf,norm3D)FwdNormMeasSum(aRecon2,ftRecons,myIllum,myOtf,norm3D,FwdModel); % normalizes the result of the FwdModel befor calculating the error and residuum
    BwdModel = @(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D) BwdNormMeasSum(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D,BwdModel); % accounts for FwdNormcomp in the calculation of the gradient
end

if RegularisationParameters(17,1)  % NormMeasSumSqr, normalizes the result of the FwdModel befor calculating the error and residuum    
    if RegularisationParameters(16,1)
        error('Only one of NormMeasSum or NormMeasSumSqr can be used.');
    end
    FwdModel=@(aRecon2,ftRecons,myIllum,myOtf,norm3D)FwdNormMeasSumSqr(aRecon2,ftRecons,myIllum,myOtf,norm3D,FwdModel); % normalizes the result of the FwdModel befor calculating the error and residuum
    BwdModel = @(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D) BwdNormMeasSumSqr(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D,BwdModel); % accounts for FwdNormcomp in the calculation of the gradient
end

if RegularisationParameters(12,1) ~= 0.0 % Background is included in . Nothing needs to be done for the Bwd direction
    FwdModel=@(aRecon2,ftRecons,myIllum,myOtf,norm3D)ApplyBgModel(aRecon2,ftRecons,myIllum,myOtf,norm3D,FwdModel,RegularisationParameters(12,1));  % fixes the second parameter is the current Fwd model and the third one to the set Background
end

%% This is a trick to speed things up for Gaussian noise model under some conditions:  Ptochography-like updates
useFTComparison=0;
if RegularisationParameters(10,1)==0  %    'LeastSqr'
    global DeconvMask
    if  strcmp(char(FwdModel),'FwdObjConvPSF') && strcmp(char(BwdModel),'BwdResidObjConvPSF')
        if isempty(DeconvMask)
            fprintf('Found the opportunity to use fast (Ptychography-like) updates. Reassigning Fwd and Bwd models.\n')
            FwdModel=@FwdHalfObjConvPSF;
            BwdModel=@BwdHalfResidObjConvPSF;
            CalcResiduum=@ResidRFTLeastSqr;
            if isempty(myFTim)
                fprintf('Transforming Data to Fourier space\n');
                for n=1:numel(myim)
                    myFTim{n}=rft(myim{n});
                end
            end
            useFTComparison=1;
        else
            fprintf('Warning: Without extra boundaries, you would be able to use fast Ptychography updates using the LeastSqr norm.\n')            
        end
    elseif  strcmp(char(FwdModel),'FwdObjIlluConvPSF') && strcmp(char(BwdModel),'BwdResidObjIlluConvPSF')
        if isempty(DeconvMask)
            fprintf('Found the opportunity to use fast (Ptychography-like) updates. Reassigning Fwd and Bwd models.\n')
            FwdModel=@FwdHalfObjIlluConvPSF;
            BwdModel=@BwdHalfResidObjIlluConvPSF;
            CalcResiduum=@ResidRFTLeastSqr;
            if isempty(myFTim)
                fprintf('Transforming Data to Fourier space\n');
                for n=1:numel(myim)
                    myFTim{n}=rft(myim{n});
                end
            end
            useFTComparison=1;
        else
            fprintf('Warning: Without extra boundaries, you would be able to use fast Ptychography updates using the LeastSqr norm.\n')            
        end
    end
end
