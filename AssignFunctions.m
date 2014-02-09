% AssignFunctions(RegulisierungsParameters,ToEstimate) : Assigns the helper functions FwdModel, BwdModel and CalcResiduum, depending on the model desctibed by ToEstimate
% ToEstimate = 0 : Object only estimation
% ToEstimate = 1 : Object estimation in a model including illumination
% ToEstimate = 2 : Illumination estimation
% ToEstimate = 3 : PSF estimation
function AssignFunctions(RegularisationParameters,ToEstimate)

global FwdModel; % A function (pointer) which is assigned outside. Options are FwdObjConvPSF() FwdObjConvASF() FwdObjIlluConfPSF() and FwdObjIlluConfASF()
global BwdModel;  % This performs the convolution of the residuum with the psf. Options are BwdModel() are BwdResidObjConvPSF() BwdResidObjConvASF() BwdResidObjIlluConvPSF() BwdResidObjIlluConvASF() BwdResidIlluObjConvPSF() BwdResidIlluObjConvASF()
global CalcResiduum; % A function (pointer) which is assigned outside. Options are ResidPoisson(), ResidLeastSrq(), ResidWeightedLeastSqr()
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF
global ConvertGradToVec; % Contains the routine to convert the model (e.g. the gradient of the object) into the vector to iterate
global ConvertModelToVec; % Contains the routine to convert the model (e.g. the gradient of the object) into the vector to iterate
global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
global myim; % This is the measured data. It is needed to see, if it is of type complex, in which case even real reconstructions have to generate full complex data.
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
else
    error('Unknown residuum calculation method');
end
FwdModel={};BwdModel={};ConvertInputToModel={};ConvertModelToVec={};ConvertGradToVec={};
my_sumcond={};
if ToEstimate==0  % Object is estimated, illumination and psf are assumed to be known
    AssignToGlobal=@AssignToObject;
    if ~isreal(myim{1}) || ComplexPSF % data is of type complex or the "IntensityData" flag is activated
        FwdModel=@FwdObjConvASF;
        BwdModel=@BwdResidObjConfASF;
        if RegularisationParameters(7,1) % case 'Complex'
            ConvertInputToModel=@convertVecToCpxObj;  % This is a packed complex vector
            ConvertModelToVec=@convertCpxObjToVec;
        else
            ConvertInputToModel=@convertRVecToCpxObj;  % This is only a real valued vector but expanded to complex.
            ConvertModelToVec=@convertCpxObjToRVec;
        end
    else  % data is real valued. Thus PSF is also real
        if RegularisationParameters(7,1) % 'Complex' reconstruction but real valued data
            if RegularisationParameters(8,1) % 'Intensity' data. In this case a complex reconstruction can be attempted
                FwdModel=@FwdObjConvASFSqr;
                BwdModel=@BwdResidObjConfASFSqr;
                ConvertInputToModel=@convertVecToCpxObj;  % This is a packed complex vector
                ConvertModelToVec=@convertCpxObjToVec;
            else
                error('Flag Complex does not make sense for real valued data, if the Intensity flag is not set.');
            end
        else % real reconstruction and real data
            ConvertInputToModel=@convertVecToObj;
            FwdModel=@FwdObjConvPSF;
            BwdModel=@BwdResidObjConfPSF;
            ConvertModelToVec=@convertObjToVec;
        end
    end
elseif ToEstimate==1  % Object is estimated, illumination is assumed known, but spatially variing
    AssignToGlobal=@AssignToObject;
    if RegularisationParameters(7,1) % case 'Complex'
        ConvertInputToModel=@convertVecToCpxObj;
        FwdModel=@FwdObjIlluConvASF;
        BwdModel=@BwdResidObjIlluConfASF;
        ConvertModelToVec=@convertCpxObjToVec;
    else
        ConvertInputToModel=@convertVecToObj;
        FwdModel=@FwdObjIlluConvPSF;
        BwdModel=@BwdResidObjIlluConfPSF;
        ConvertModelToVec=@convertObjToVec;
    end
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
        BwdModel=@BwdResidIlluObjConfASF;
    else
        FwdModel=@FwdObjIlluConvPSF;
        BwdModel=@BwdResidIlluObjConfPSF;
    end
elseif ToEstimate==3  % object and illumination are assumed to be known, and otf is estimated
    if RegularisationParameters(14,1) % case 'ProjPupil' where only the 2D pupil is estimated
        AssignToGlobal=@AssignToOTF;
        ConvertInputToModel=@convertPupilVecToOTF;
        FwdModel=@FwdObjConvPSF;   % identical for object and PSF (already converted to OTF)
        BwdModel=@BwdResidPSFConfObj;  % has to convolve with aRecon
        ConvertModelToVec=@convertPSFToPupilVec;  % this conversion is identical for object or PSF
    else
        AssignToGlobal=@AssignToOTF;
        ConvertInputToModel=@convertVecToOTF;  % fills a volume in Fourier-space
        FwdModel=@FwdObjConvPSF;   % identical for object and PSF (already converted to OTF)
        BwdModel=@BwdResidPSFConfObj;  % has to convolve with aRecon
        ConvertModelToVec=@convertPSFToOTFVec;  % this conversion is identical for object or PSF
    end
end

if isempty(ConvertGradToVec);
    ConvertGradToVec=ConvertModelToVec;    
end

if RegularisationParameters(9,1) % ForcePos, Forces object to be positive by estimating only its squareroot
    ConvertInputToModel=@(vec)ApplySqrModel(vec,ConvertInputToModel);  % fixes the second parameter is the current model
    ConvertGradToVec=@(grad)ApplySqrModelGrad(grad,ConvertGradToVec);
    ConvertModelToVec=@(model)ApplyInverseSqrModel(model,ConvertModelToVec);
end

if RegularisationParameters(12,1) ~= 0.0 % Background is included in . Nothing needs to be done for the Bwd direction
    FwdModel=@(aRecon2,ftRecons,myIllum,myOtf,norm3D)ApplyBgModel(aRecon2,ftRecons,myIllum,myOtf,norm3D,FwdModel,RegularisationParameters(12,1));  % fixes the second parameter is the current Fwd model and the third one to the set Background
end
