% [err,thegrad]=GenericErrorAndDeriv(myinput) : Error measure for deconvolution algorithm

%***************************************************************************
%   Copyright (C) 2008-2009 by Rainer Heintzmann                          *
%   heintzmann@gmail.com                                                  *
%                                                                         *
%   This program is free software; you can redistribute it and/or modify  *
%   it under the terms of the GNU General Public License as published by  *
%   the Free Software Foundation; Version 2 of the License.               *
%                                                                         *
%   This program is distributed in the hope that it will be useful,       *
%   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
%   GNU General Public License for more details.                          *
%                                                                         *
%   You should have received a copy of the GNU General Public License     *
%   along with this program; if not, write to the                         *
%   Free Software Foundation, Inc.,                                       *
%   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
%**************************************************************************
%

function [err,thegrad]=GenericErrorAndDeriv(myinput)
global myim;
global myFTim;
global myillu;  % only of there is an illumination pattern present, will it be used
global otfrep;
global BwdOTF;
global BetaVals;  % These are the scaling factors between pixel coordinates. This is proportional to the pixel width, but should be normalized
global DeconvMask;  % only data in this mask will be evaluated
global ToEstimate;   % This flag controls what to estimate in this iteration step. 0: sample density, 1: illumination intensity, 2: both, 3: psf
global aRecon;   % Here the sample is stored, if the estimation is illumination or psf.
global aResampling;   % If unequalt to one, a resampling will be introduced betwenn reconstruction and measured data.
global subSampling;   % If unequalt to one, a resampling will be introduced betwenn reconstruction and measured data.
global NormFac;  % normalisation factor
global my_sumcond;   % denotes the positions in myillu for which each the sum condition (down to previous mention) are fullfilled.

global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF
global ConvertGradToVec; % Contains the routine to convert the model (e.g. the gradient of the object) into the vector to iterate
global CalcResiduum; % A function (pointer) which is assigned outside. Options are ResidPoisson(), ResidLeastSrq(), ResidWeightedLeastSqr()
global FwdModel; % A function (pointer) which is assigned outside. Options are FwdObjConvPSF() FwdObjConvASF() FwdObjIlluConfPSF() and FwdObjIlluConfASF()
global BwdModel;  % This performs the convolution of the residuum with the psf. Options are BwdModel() are BwdResidObjConvPSF() BwdResidObjConvASF() BwdResidObjIlluConvPSF() BwdResidObjIlluConvASF() BwdResidIlluObjConvPSF() BwdResidIlluObjConvASF()
global BwdExtraModel;  % Here changes to the Bwd projection are applied to, such as Sqr or Phase
global RegularisationParameters;  % Used only for 'Show'
global measSumsSqr;
global measSums;
global measSum;
global measSumSqr;
global useFTComparison;  % only active if Ptychography updates are used in this round
global Recons;  % to make it accessible from the outside
global RegObjVal; % stores the object regularization, such that it is also available as a constant for the energy term during PSF estimation
global savedATF; % A java viewer in which to show the Deconv progress continuously.

global DeconvViewer; % A java viewer in which to show the Deconv progress continuously.
global IllumViewer; % A java viewer in which to show the Illumination estimation progress continuously.
global PSFViewer; % A java viewer in which to show the PSF progress continuously.
global OTFViewer; % A java viewer in which to show the OTF progress continuously.
global ATFViewer; % A java viewer in which to show the pupil progress continuously.

%delta= 100; % Weight for the negativtiy penalty
%delta= 1000; % Weight for the negativtiy penalty
%delta= 1; % Weight for the negativtiy penalty

DMask=[];
ftRecon=[];
if ~isempty(DeconvMask)
    DMask=DeconvMask;
end

% Recast the matlab data into the dip_image datastructure
%aRecon=dip_image(aRecon);
norm3D = sqrt(prod(to3dvec(size(myim{1})) .* to3dvec(subSampling)));  % aRecon is not defined yet.
norm3DObj = sqrt(prod(floor(to3dvec(aResampling) .* to3dvec(subSampling) .* to3dvec(size(myim{1})))));  % aRecon is not defined yet.    .* to3dvec(subSampling)
% DataSize=size(myim{1});
% DataLength=prod(DataSize);
%if length(DataSize) < 3
%    DataSize(3) = 1;
%end
%% First fill in Object estimate, illumination or PSF from the Matlab type input vector

thegrad=AssignToGlobal(ConvertInputToModel(myinput)); % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF
% However, it can also be several functions wrapped together to allow for taking the absolute square (to force positivity)
% ConvertInputToModel is a global function, assigned to a specific function
% It returns an empty gradient vector (if asked to do so)
% It also accounts for the ForcePositiv constraint

if ToEstimate==2  % PSF/OTF estimation
    thegrad=repmat(thegrad,[1 1 1 numel(otfrep)]);
end

err=0;           % clears the errorsum
clear myinput;

% dip_setboundary('periodic')  % Establishes Periodic Boudary conditions.

currentSumCondIdx=1;   % goes through each of the sum conditions
prevSumCondGradIdx=0;

if isempty(ToEstimate) || ToEstimate==0  % Object estimate is only regularised once even for multi-view deconv
    [myReg,myRegGrad]=Regularize(aRecon,BetaVals);  % Object regularisation is applied only once!
    RegObjVal=myReg;  % stores this for potential use during PSF optimization. WIHTOUT NORMFAC
    err=err+myReg;    % was cleared before the for-loop
    % myRegGrad=BwdExtraModel(myRegGrad,1);  % The background model is applied after the for loop for all instances    
    thegrad=thegrad+ myRegGrad;
    if RegularisationParameters(13,1)
        DeconvViewer=LiveView(DeconvViewer,aRecon);
    end
else
    if ~isempty(RegObjVal)
        err=err+RegObjVal;  % for potential use during PSF optimization. In this case there is no gradient! WITHOUT NormFac!
    end
end

ftRecon=[];
numViews=length(myim);
OTFsPerView=floor(length(myillu)/length(myim));  % e.g. one illu per image means OTFsPerView=1
if OTFsPerView < 1
    OTFsPerView = 1;
end
% if OTFsPerView ~= size(myillu,1)
if OTFsPerView ~= size(myillu,1) && (~isempty(myillu)) %Aurelie 03.06.2014 for the WF case
    error('illumination patterns dont agree with number of OTFs')
end

for viewNum = 1:numViews    % This iterates over all the different measured images
    myReg=0;myRegGrad=0;

    
    %% select the appropriate subdata regions
    if ~useFTComparison
        myImg=myim{viewNum};  % This does not cost any time or memory
    else
        myImg=myFTim{viewNum};  % The Ptychography trick is used to avoid some FTs. Thus the data has to be compared with the Fourier transforms
    end
    measSum=measSums{viewNum};
    measSumSqr=measSumsSqr{viewNum};

    for subViewOTFNum =1:OTFsPerView  % iterates over sup-views in the case of 3DSIM generating as a sum only one Fwd projected image
        
        myOtf=otfrep{1+mod(viewNum-1+(subViewOTFNum-1),length(otfrep))};  % This does not cost any time or memory. If only one otf is there it will always be used
        if ~isempty(myillu)
%             myIllum=myillu{(subViewOTFNum-1),1+mod((viewNum-1),length(myillu))};
            myIllum=myillu{subViewOTFNum,1+mod((viewNum-1),length(myillu))}; %Aurelie 26.05.2014
        else
            myIllum=1;
        end
        
        if viewNum==1 && (isempty(ToEstimate) || ToEstimate==0 || ToEstimate==2) && (isempty(myillu))% Object or OTF estimate
            if ~isequal(FwdModel,@FwdIdentity)  % Otherwise the ft is not needed.
                if isreal(aRecon) && ~equalsizes(size(myOtf),size(aRecon))
                    ftRecon=rft(aRecon);  % To save time and only compute this ft once for multi-view deconvolutions
                    if any(aResampling~=1)
                        ftRecon=rft_resize(ftRecon,1./aResampling);  % if the user wants to use a different reconstruction grid
                    end
                else
                    ftRecon=ft(aRecon);  % To save time and only compute this ft once for multi-view deconvolutions
                end
            end
        end
        %% first we need to apply the forward model
        if subViewOTFNum == 1
            % also serves as an initialization for the sum
            Recons = FwdModel(aRecon,ftRecon,myIllum,myOtf,norm3D^2/norm3DObj);  % This performs the convolution of the object (multiplied with illumination) with the psf
        else
            Recons = Recons+FwdModel(aRecon,ftRecon,myIllum,myOtf,norm3D^2/norm3DObj);  % This performs the convolution of the object (multiplied with illumination) with the psf
        end
    end  % of iteration over sub-views (needed for 3D SIM forward model)
    Recons=FwdSubsample(Recons,subSampling);
    % Functions to be hidden behind FwdModel() are FwdObjConvPSF() FwdObjConvASF() FwdObjIlluConfPSF() and FwdObjIlluConfASF()
    %% Now calculate the residuum depending on the deconvolution method
    [myError,residuum] = CalcResiduum(Recons,myImg,DMask); % May change Recons to force positivity. Possible Functions are ResidPoisson, ResidLeastSqr, ResidWeightedLeastSqr
    % clear Recons;
    err=err+myError;
    residuum=BwdSubsample(residuum,subSampling);
    %% Apply the Transpose (Bwd Model)
    for subViewOTFNum =1:OTFsPerView  % iterates over sup-views in the case of 3DSIM generating as a sum only one Fwd projected image
        OTFNum = 1+mod(viewNum-1+(subViewOTFNum-1),length(otfrep));
        myOtf=otfrep{OTFNum};  % This does not cost any time or memory. If only one otf is there it will always be used
        if ~isempty(myillu)
%             myIllum=myillu{(subViewOTFNum-1),1+mod((viewNum-1),length(myillu))};
            myIllum=myillu{subViewOTFNum,1+mod((viewNum-1),length(myillu))}; %Aurelie 26.05.2014
        else
            myIllum=1;
        end

        if RegularisationParameters(35,1) && ToEstimate~=2
            BwdOtf = BwdOTF{OTFNum};
        else
            BwdOtf = myOtf;
        end
        
        if subViewOTFNum == 1
            % also serves as an initialization for the sum
            myGrad=BwdModel(residuum,aRecon,ftRecon,myIllum,BwdOtf,norm3DObj);  % This performs the convolution of the residuum with the psf
        else
            myGrad=myGrad+BwdModel(residuum,aRecon,ftRecon,myIllum,BwdOtf,norm3DObj);  % This performs the convolution of the residuum with the psf
        end
        
        if viewNum==numViews
            clear ftRecon;
        end
        if ToEstimate ~= 2
            clear myOtf;
        end
        clear residuum;
        % Functions to be hidden behind BwdModel() are BwdResidObjConvPSF() BwdResidObjConvASF() BwdResidObjIlluConvPSF() BwdResidObjIlluConvASF() BwdResidIlluObjConvPSF() BwdResidIlluObjConvASF()
        
        if isempty(ToEstimate) || ToEstimate==0  % object estimate has to be summed for all views
            thegrad=thegrad + myGrad;
            % myGrad=BwdExtraModel(myGrad,1);  % The background model is applied after the for loop for all instances
            % thegrad=thegrad + myGrad;
            % clear myGrad;
        else
            myGrad=BwdExtraModel(myGrad,viewNum);  % Possibly changes due to Sqr or Phase  ('ForcePos' or 'ForcePhase')
        end
        clear myImg;
        
    end % of iteration over sub-views (needed for 3D SIM forward model).. thegrad is summed in this loop
    if isempty(ToEstimate) || ToEstimate==0  % object estimate has to be summed for all views
        clear myGrad;
    end
    % The second part of the gradient calculation (multiplication with aRecon or illu{viewNum} is done at the end of the loop    
    % The terms below are only dependend on the reconstruction but not on the data. They need to be calculated onyl once for all views.
    % Total error term (including regularisation):
    if ~isempty(ToEstimate) && ToEstimate>0  % illumination or OTF estimate, once for each view
        if ToEstimate==1
            [myReg,myRegGrad]=Regularize(myIllum,BetaVals);
            if RegularisationParameters(13,1)
                IllumViewer=LiveView(IllumViewer,myIllum);
            end
        else
            if NeedsRegularisation()
                mypsf=ifftshift(rift(myOtf));
                %mypsf=circshift(mypsf,floor(size(mypsf)/2));
                [myReg,myRegGrad]=Regularize(mypsf,BetaVals);  % Should this better be the PSF ?
                % myRegGrad=circshift(myRegGrad,-floor(size(mypsf)/2));  % back to the coordinate system of the gradient
                myRegGrad=fftshift(myRegGrad);
                clear mypsf;
            end
        end
        agradIdx=viewNum-1-(currentSumCondIdx-1);
        if ToEstimate==2 || isempty(my_sumcond) || viewNum ~= my_sumcond{currentSumCondIdx}  % OTF estimation or illumination estimation without a sumcondition
            if ToEstimate==2 && (agradIdx >= length(otfrep))
                error('Number of images does not correspond to number of OTFs for blind OTF deconvolution');
            end
            thegrad(:,:,:,agradIdx)=myGrad + myRegGrad;  % thegrad estimates only ONE pattern even if there are multiple sub-patterns in it
        else  % The last residuum has to be subtracted from each of the other residuals, see eq. S14 and S4 in supplementary methods of DOI: 10.1038/NPHOTON.2012.83
            thegrad(:,:,:,prevSumCondGradIdx:agradIdx-1)=thegrad(:,:,:,prevSumCondGradIdx:agradIdx-1)-(myGrad+myRegGrad);
            % thegrad(:,:,:,prevSumCondGradIdx:agradIdx-1)=thegrad(:,:,:,prevSumCondGradIdx:agradIdx-1)-repmat(myGrad+myRegGrad,[1 1 1 agradIdx-prevSumCondGradIdx]);
            prevSumCondGradIdx=agradIdx;
            currentSumCondIdx = currentSumCondIdx +1;
        end
        clear myGrad;
        
        err=err+myReg;
    end

end  % viewNum
if ToEstimate == 2 && RegularisationParameters(13,1)
    allotf=cat(4,otfrep{:});
    PSFViewer=LiveView(PSFViewer,ifftshift(rift3d(allotf)));
    OTFViewer=LiveView(OTFViewer,allotf);
    if ~isempty(savedATF)
        ATFViewer=LiveView(ATFViewer,cat(1,savedATF(:,:,:,0),savedATF(:,:,:,1),savedATF(:,:,:,2)));
    end
end

if isempty(ToEstimate) || ToEstimate==0  % object estimate has to be summed for all views
   thegrad = BwdExtraModel(thegrad,1);
end
clear DMask;

err = double(NormFac*err);

thegrad = NormFac * ConvertGradToVec(thegrad);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation

