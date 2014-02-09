% MyIdivErrorAndDeriv(aRecon) : Error measure for deconvolution algorithm
% this function interpretes the 4th dimension as a multi-view deconvolution
% each element of the 4th dimension having a corresponding PSF

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
global myillu;  % only of there is an illumination pattern present, will it be used
global otfrep;
global BetaVals;  % These are the scaling factors between pixel coordinates. This is proportional to the pixel width, but should be normalized
global DeconvMask;  % only data in this mask will be evaluated
global ToEstimate;   % This flag controls what to estimate in this iteration step. 0: sample density, 1: illumination intensity, 2: both, 3: psf
global aRecon;   % Here the sample is stored, if the estimation is illumination or psf.
global aResampling;   % If unequalt to one, a resampling will be introduced betwenn reconstruction and measured data.
global NormFac;  % normalisation factor
global my_sumcond;   % denotes the positions in myillu for which each the sum condition (down to previous mention) are fullfilled.

global AssignToGlobal; % Assigns the converted data to the appropriate global variable and return an empty gradient vector
global ConvertInputToModel; % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF
global ConvertGradToVec; % Contains the routine to convert the model (e.g. the gradient of the object) into the vector to iterate
global CalcResiduum; % A function (pointer) which is assigned outside. Options are ResidPoisson(), ResidLeastSrq(), ResidWeightedLeastSqr()
global FwdModel; % A function (pointer) which is assigned outside. Options are FwdObjConvPSF() FwdObjConvASF() FwdObjIlluConfPSF() and FwdObjIlluConfASF()
global BwdModel;  % This performs the convolution of the residuum with the psf. Options are BwdModel() are BwdResidObjConvPSF() BwdResidObjConvASF() BwdResidObjIlluConvPSF() BwdResidObjIlluConvASF() BwdResidIlluObjConvPSF() BwdResidIlluObjConvASF()
global RegularisationParameters;  % Used only for 'Show'

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
norm3D = sqrt(prod(size(myim{1})));  % aRecon is not defined yet.
DataSize=size(myim{1});
DataLength=prod(DataSize);
%if length(DataSize) < 3
%    DataSize(3) = 1;
%end
%% First fill in Object estimate, illumination or PSF from the Matlab type input vector

thegrad=AssignToGlobal(ConvertInputToModel(myinput)); % Will change either aRecon, myillu or otfrep. It can be convertVecToObj, convertVecToIllu, convertVecToPSF
% However, it can also be several functions wrapped together to allow for taking the absolute square (to force positivity)
% ConvertInputToModel is a global function, assigned to a specific function
% It returns an empty gradient vector (if asked to do so)
% It also accounts for the ForcePositiv constraint

if ToEstimate==2
    thegrad=repmat(thegrad,[1 1 1 numel(otfrep)]);
end

err=0;           % clears the errorsum
clear myinput;

% dip_setboundary('periodic')  % Establishes Periodic Boudary conditions.

currentSumCondIdx=1;   % goes through each of the sum conditions
prevSumCondGradIdx=0;

if isempty(ToEstimate) || ToEstimate==0  % Object estimate is only regularised once even for multi-view deconv
    [myReg,myRegGrad]=Regularize(aRecon,BetaVals);  % Object regularisation is applied only once!
    err=err+myReg;    % was cleared before the for-loop
    thegrad=thegrad+ myRegGrad;
    if RegularisationParameters(13,1)
        dipshow(3,aRecon);drawnow();
    end
end

ftRecon=[];
numViews=length(myim);
for viewNum = 1:numViews    % This iterates over all the different measured images
    myReg=0;myRegGrad=0;

    
    %% select the appropriate subdata regions
    myImg=myim{viewNum};  % This does not cost any time or memory
    myOtf=otfrep{1+mod(viewNum-1,length(otfrep))};  % This does not cost any time or memory. If only one otf is there it will always be used
    if ~isempty(myillu)
        myIllum=myillu{1+mod(viewNum-1,length(myillu))};
    else
        myIllum=1;
    end
    
    if viewNum==1 && (isempty(ToEstimate) || ToEstimate==0 || ToEstimate==2) % Object or OTF estimate
        if isreal(aRecon) && (norm(size(myOtf)-size(aRecon))~=0)
            ftRecon=rft(aRecon);  % To save time and only compute this ft once for multi-view deconvolutions
            if any(aResampling~=1)
                ftRecon=rft_resize(ftRecon,1./aResampling);  % if the user wants to use a different reconstruction grid
            end
        else
            ftRecon=ft(aRecon);  % To save time and only compute this ft once for multi-view deconvolutions
        end
    end
    %% first we need to apply the forward model
    Recons = FwdModel(aRecon,ftRecon,myIllum,myOtf,norm3D);  % This performs the convolution of the object (multiplied with illumination) with the psf
    % Functions to be hidden behind FwdModel() are FwdObjConvPSF() FwdObjConvASF() FwdObjIlluConfPSF() and FwdObjIlluConfASF()
    %% Now calculate the residuum depending on the deconvolution method
    [myError,residuum] = CalcResiduum(Recons,myImg,DMask); % May change Recons to force positivity. Possible Functions are ResidPoisson, ResidLeastSqr, ResidWeightedLeastSqr
    clear Recons;
    err=err+myError;
    %% Apply the Transpose (Bwd Model)
    myGrad=BwdModel(residuum,aRecon,ftRecon,myIllum,myOtf,norm3D);  % This performs the convolution of the residuum with the psf
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
        clear myGrad;
    end
    clear myImg;
    
    % The second part of the gradient calculation (multiplication with aRecon or illu{viewNum} is done at the end of the loop    
    % The terms below are only dependend on the reconstruction but not on the data. They need to be calculated onyl once for all views.
    % Total error term (including regularisation):
    if ~isempty(ToEstimate) && ToEstimate>0  % illumination or OTF estimate, once for each view
        if RegularisationParameters(13,1)
            dipshow(3,myIllum);drawnow();
        end
        if ToEstimate==1
            [myReg,myRegGrad]=Regularize(myIllum,BetaVals);
        else
            mypsf=ifftshift(rift(myOtf));
            %mypsf=circshift(mypsf,floor(size(mypsf)/2));
            [myReg,myRegGrad]=Regularize(mypsf,BetaVals);  % Should this better be the PSF ?
            % myRegGrad=circshift(myRegGrad,-floor(size(mypsf)/2));  % back to the coordinate system of the gradient
            myRegGrad=fftshift(myRegGrad);
            clear mypsf;
        end
        agradIdx=viewNum-1-(currentSumCondIdx-1);
        if ToEstimate==2 || isempty(my_sumcond) || viewNum ~= my_sumcond{currentSumCondIdx}  % OTF estimation or illumination estimation without a sumcondition
            if ToEstimate==2 && (agradIdx >= length(otfrep))
                error('Number of images does not correspond to number of OTFs for blind OTF deconvolution');
            end
            thegrad(:,:,:,agradIdx)=myGrad + myRegGrad;
        else  % The last residuum has to be subtracted from each of the other residuals, see eq. S14 and S4 in supplementary methods of DOI: 10.1038/NPHOTON.2012.83
            thegrad(:,:,:,prevSumCondGradIdx:agradIdx-1)=thegrad(:,:,:,prevSumCondGradIdx:agradIdx-1)-repmat(myGrad+myRegGrad,[1 1 1 agradIdx-prevSumCondGradIdx]);
            prevSumCondGradIdx=agradIdx;
            currentSumCondIdx = currentSumCondIdx +1;
        end
        clear myGrad;
        
        err=err+myReg;
    end
    
end  % viewNum
clear DMask;

err = double(NormFac*err);

thegrad=ConvertGradToVec(NormFac*thegrad);    % converts the dip_image back to a linear matlab vector. Also does the required Fourier-transform for illumination estimation

