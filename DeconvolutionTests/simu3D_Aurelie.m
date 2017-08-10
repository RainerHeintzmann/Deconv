%% Multi-view blind SIM deconvolution
% Create an object that demonstrates the optical sectioning power as well as the resolution of blind-SIM
useCuda=0; disableCuda();
addpath('C:\Users\user\Documents\MATLAB\Aurelie')


OldVersion=0; %If 1, uses the back-up version
withAbberation=0; %If 1, uses an abberated illumination
dbg_show=0; %Shows a number of intermediate results
% switchVersion(OldVersion) % in=0: version on the server. in=1: back-up from 12.09.2013
if ~OldVersion
    SuperSampling=2;
end
%% Create object
% a=readim('Y:\PublicFiles\MATLAB\images\resolution_512.tif');
% a=readim('Y:\MATLAB\images\resolution_fine.tif');
% sV=[256 256]; %Goal size
% convX=1/floor(size(a,1)./sV(1)); %conversion of size (zoom to apply) in X
% convY=1/floor(size(a,2)./sV(2)); %conversion of size (zoom to apply) in Y
% % a2=resample(a,[0.5 0.5],[0 0]); %Resample image (factor 2)
% a2=resample(a,[convX convY],[0 0],'linear'); %Resample image. Ugly

% b=rotation(a,pi);
% c=rotation(a,pi/2);

if (1) %spoke with out-of-focus bar
    sX=200;sY=200;sZ=23;
    NumSpokes=40; %For the Spoke object. 40
    zFoc=3; % Slice index of the in-focus information
    a=GenObj('obj','Spoke','sv',[sX sY sZ], 'NumSpokes', NumSpokes,'Zfoc',zFoc);
    obj=a; %rename to avoid confusion later
else
    a=readim('Y:\MATLAB\images\resolution_512.tif'); %whole resolution target
    d=rotation(a,3*pi/2);
    if OldVersion
        obj=newim([size(a,1),size(a,2),8]); %Old version
    else
    %     obj=newim([size(a,1),size(a,2),9]); %New
    %     obj=newim([size(a,1),size(a,2),23]); %Rainer, for good sampling
        obj=newim([size(a,1),size(a,2),13]); %To stop get confused about nb of images and nb of planes
    end
    zFoc=3; %focal slice
    % dist=2; %distance between 2 objects
    obj(:,:,zFoc)=a; %obj(:,:,zFoc+1)=b*1.5;
    %obj(:,:,zFoc+2)=c*1.5; 
    % obj(:,:,zFoc+3)=d*1.5; %600 nm
    obj(:,:,zFoc+4)=d*1.5; %800 nm distance
    % obj(:,:,zFoc+dist)=b*1.5;
    clear b
    clear c
    clear d
    a=obj; %rename to avoid confusion later
end



if dbg_show
    a
end

%% PSF
% For the multi-view (3D) case, we have to simulate 2 OTFs
% TO DO: see how it is done in Simon's code
% TO DO: figure out how to pass it as an argument through the input which is a PSF

 useCuda=0; disableCuda();

% Parameters

sX=size(obj,1); sY=size(obj,2); sZ=size(obj,3);
scaleX=80; scaleY=scaleX; scaleZ=100; % pixel sizes in nm (along z it is always controversary)
NA=1.3;
lambda=500;
ri=1.52;

% SIMDist=lambda/(1.4*NA)/scaleX;  % Emerics Simulation
SIMDist=2*lambda/(1.8*NA)/scaleX; % at 90% to cutoff for fine grating. Coarse is given here.
% SIMDist=lambda/(2.8*NA)/scaleX;  % higher Frequency (beyond cutoff)
NumPhotons=100000 * (4*NA)^2/scaleX/scaleY;    % Photon number in Nyquist pixels
% NumPhotons=100 * (4*NA)^2/scaleX/scaleY;    % Photon number in Nyquist pixels
NumPhotons=NumPhotons*100;
Offset=0;


if (1)
    ImageParam=struct('Sampling',[scaleX scaleY scaleZ],'MaxPhotons',0,'Size',[size(obj)]); %Not sure MaxPhotons is needed
    PSFParam=struct('NA',NA,'n',ri,'MinOtf',1.2e-3,'lambdaEm',lambda); %Not sure MinOtf is needed
    k0=[sX 0]/(SIMDist); %this should be correct
    OTFRelPhase=[0 0]; %should have 2 components if 2 OTFs are needed
%     otfs=GenOTFs(ImageParam,PSFParam,k0,OTFRelPhase,minotf)
    otfs=GenOTFs(ImageParam,PSFParam,k0,OTFRelPhase); %Use the default minotf value 1.2e-3
    h=cell(numel(otfs),1);
    for v=1:numel(otfs)
        h{v}=fftshift(rift(otfs{v}));
    end
    %now, h is a cell array of 2 components with 2 PSFs
else
    h=kSimPSF({'sX',sX;'sY',sY;'sZ',sZ;'scaleX',scaleX;'scaleY',scaleY;'scaleZ',scaleZ;'confocal',0;'na',NA;'lambdaEm',lambda;'ri',ri});
    % writeim(h,'psf_3D_scaleZ50') %30 slices: we can cut later
    % h=readim('psf_3D_scaleZ50');
    % h=extract(h,[sX sY 20]); %20 slices are enough
    h2=fixWFPSF(h,0.4); %the default value for the parameter is 0,35. Here, not enough.
    % Remark: the bigger the size, the higher this parameter should be taken
    % Or smaller scaleZ
    % h2=fixWFPSF(h,0.99); %always check, sometimes
    sum(h,[],[1 2]); %Is not flat (i.e., the sampling of the PSF leads to higher total intensity in the focal plane)
     sum(h2,[],[1 2]); %corrected. Now, the WF img of a plane is not resolvable anymore
     if(1)
         h=h2; %use the correct PSF
     end
     if dbg_show
         sum(h,[],[1 2])
     end
end



%% SIM images

% mean(h,[],[1 2]); %should be flat
% real(ift(ft(obj).*ft(h))); %check that in the X-Z view, the flattening is visible
% idwf=real(ift(ft(obj).*ft(h2))); %check that in the X-Z view, the flattening is visible
% idwf(:,:,zFoc);

% AbberationMap= (xx(obj)*xx(obj))/(100*100) *pi;
if withAbberation
    map=abs(xx(size(a))) < 70;
    AbberationMap=-gaussf(map,3)*pi; %That's the one!
    myGrating=cell(numel(otfs),9);
%     [myillu,myamp]=GenSIM(asize,D,numphases,numdirs,contrast,startdir, AbberationMap)
    myGrating{1,:}=GenSIM(size(a),SIMDist,5,3,1,AbberationMap); %Illu pattern
    myGrating{2,:}=GenSIM(size(a),SIMDist/2,5,3,1,AbberationMap); %SIMDist finest grating period???
    myGratingTh=cell(numel(otfs),9);
    myGratingTh{1,:}=GenSIM(size(a),SIMDist,5,3,1,0); %Illu pattern non distorted
    myGratingTh{2,:}=GenSIM(size(a),SIMDist/2,5,3,1,0); %SIMDist finest grating period???
%     myGratingTh=GenSIM(size(a),SIMDist,3,3,1,0); %Illu pattern non distorted
else
    AbberationMap=0;
%     myGrating=cell(numel(otfs),9);
%      myGrating=cell(numel(otfs),1);
%     [myillu,myamp]=GenSIM(asize,D,numphases,numdirs,contrast,startdir, AbberationMap)
    aGrating1=GenSIM(size(a),SIMDist,5,3,1,AbberationMap); %Illu pattern
    aGrating2=GenSIM(size(a),SIMDist/2,5,3,1,AbberationMap); %SIMDist finest grating period???
    myGrating=[aGrating1; aGrating2]; %Combination of 2 cell arrays: myGrating is a 2*9 cell array
%     myGrating=GenSIM(size(a),SIMDist,3,3,1,AbberationMap); %Illu pattern
    myGratingTh=myGrating; %No need to recall GenSIM
    
end
% AbberationMap=-gaussf(abs(xx(size(a))) < 30,8) *pi; %Not enough aberrated?
% AbberationMap=-gaussf(abs(xx(size(a))) < 50) *pi;
% myGrating=GenSIM(size(a),SIMDist,3,3,1,0)


global myillu;

if (0) %2 beam SIM: wrong here
    % myillu=myGrating;
    myillu=aGrating1; %2 beam SIM
    %This is wrong.

    % otfs=cell(1,1); %would overwrite
    img=cell(numel(myillu),1);
    % otfs{1} = ft(h); %would overwrite: name diffently
    ImOtfs{1} = ft(h{1});
    normFac=1;
    if NumPhotons > 0   % calcualte normalization only once as speckles may change
        normFac=NumPhotons;
    end
    obj=Offset+normFac*obj/max(obj);

    for p=1:numel(myillu)
        fobj = ft(obj.*myillu{p});
        mcconv=sqrt(prod(size(obj))) * real(ift(fobj .* ImOtfs{1}));
        mcconv(mcconv<0)=0;
        if NumPhotons > 0
            img{p}=noise(mcconv,'poisson');  % put some noise on the image
        else
            img{p}=mcconv;  % put no noise
        end
    end

else
    normFac=1;
        if NumPhotons > 0   % calcualte normalization only once as speckles may change
            normFac=NumPhotons;
        end
        obj=Offset+normFac*obj/max(obj);
        
    img=cell(numel(aGrating1),1);
        
    for n=1:size(otfs,2)
        % img=img+rift2d(rft2d(repmat(obj,[1 1 1 size(illu{n},4)]) .* illu{n}) .* repmat(otfs{n},[1 1 1 size(illu,4)])); % convolve the object after multiplication with illumination with the otf
    %     myillu=repmat(myGrating{n},[1 1 size(obj,3) 1]);  % illu contains the illuminations mixed for the orders.
%         myillu=myGrating{n};
        myillu=cat(4,myGrating{n,:});
%         myotf=repmat(ft(h{n}),[1 1 1 size(myGrating{1},4)]);
%         myotf=ft(h{n}); %should be rft: this is otfs
        for p=1:size(myillu,4)
            mcconv=rift3d(rft3d(obj.* squeeze(myillu(:,:,:,p-1))) .* otfs{n}); % convolve the object after multiplication with illumination with the otf
            mcconv(mcconv<0)=0;
            if NumPhotons > 0
                mcconv=noise(mcconv,'poisson');  % put some noise on the image
            else
                mcconv=mcconv;  % put no noise
            end
            if n==1
                img{p}=mcconv;
            else
                img{p}=img{p}+mcconv;
            end

        end
    end
end

acq_imgF=cat(4,img{:});
wfimg=squeeze(mean(acq_imgF,[],4)); %Wide-field image

if dbg_show
    acq_imgF
    cat(4,myGrating{:})
 end


%%  Extract only the in-focus slice from a 3D simulation

if sZ> 1
    for n= 1:length(myGrating)
        img{n}=squeeze(img{n}(:,:,zFoc)); %In-focus slice index Zfoc (slice 3; slice 6 is here the out-of-focus info)
%         myGrating{n}=squeeze(myGrating{n}(:,:,Zfoc));
        %  Extracting the in-focus slice for the grating is ONLY needed if you want to do a two-D (non-thick slice) reconstruction
    end
    % h=squeeze(h(:,:,mid));
end

if dbg_show
    acq_img=cat(4,img{:})
 end


%% WF Deconvolution
if (0)
useCuda=1; enableCuda();
myillu=[];
img2=squeeze(mean(cat(4,img{:}),[],4)); %Wide-field image
psf2D=squeeze(h(:,:,floor(size(h,3)/2)));
myDeconvWF=GenericDeconvolution(img2,psf2D,100,'LeastSqr',[],{'ForcePos',[];'GS',0.01},[1,1],[0,0],[],useCuda)
myDeconvWF2Ds=GenericDeconvolution(img2,psf2D,100,'LeastSqr',[],{'ForcePos',[];'GS',1e-2;'Resample',[SuperSampling SuperSampling]},[1,1],[0,0],[],useCuda)
% myDeconv3DWF=GenericDeconvolution(wfimg,h,100,'LeastSqr',[],{'ForcePos',[];'GS',0.005},[1,1,3],[0,0],[],useCuda)

myDeconv3DWF=GenericDeconvolution(wfimg,h,100,'LeastSqr',[],{'ForcePos',[];'GS',0.01;'Resample',[SuperSampling SuperSampling 1]},[1,1,3],[0,0],[],useCuda);
% myDeconv3DWF=GenericDeconvolution(wfimg,h,100,'LeastSqr',[],{'ForcePos',[];'GS',0.7},[1,1,10],[0,0],[],useCuda);
myDeconv3DWF(:,:,zFoc)

end

%% Thick slice blind-SIM - preparation
 useCuda=0; disableCuda();
 
global myillu_sumcond;
% myillu_sumcond={3,6,9};  % these are the ones for which the sum condition is each fulfilled
% myillu_sumcond={6,12,18};  % necessary for 3D SIM with multi-view deconv. Doesnt work for now

% sZ=20; %Consider much less planes
% myGrating2=GenSIM([size(a,1) size(a,2) 1],SIMDist,3,3,1,0); %Just one slice
% myGrating2=myGrating;
% myGrating2=myGratingTh;
% 
% global myillu_mask;   % confines the variable only to a subspace of illu
% myillu_mask=[]; 
% myillu_mask=cell(1,numel(myGrating2));
% AberrationTolerance=15; %5 for the non-distorted case
% %dip_setboundary('periodic')  % Establishes Periodic Boudary conditions.
% for v=1:numel(myGrating2)
%     MaskThresh=20;
%     %if (useCuda)
%         %myillu_mask{v}=riftshift(cuda(dilation(dip_image_force(abs(rftshift(rft(cuda(myspecklesIdeal{v}))))> MaskThresh),AberrationTolerance)));
%     %myillu_mask{v}=riftshift(dilation(dip_image_force(abs(rftshift(rft(myspecklesIdeal{v})))> MaskThresh),AberrationTolerance));
%     myillu_mask{v}=fft2rft(ifftshift(dilation(abs(ft(myGrating2{v})) > MaskThresh,AberrationTolerance))); % Rainer's version 21.03.14
%     % myillu_mask{v}=riftshift(dilation(dip_image(abs(rftshift(rft(myGrating2{v})))> MaskThresh),AberrationTolerance));
%     %myillu_mask{v}(:,:,1:end)=0;
%     %else
%     %    myillu_mask{v}=dilation(abs(ft(myspecklesIdeal{v}))> MaskThresh,AberrationTolerance);
%     % end
% end
if withAbberation
    GenMaskFromIllu(myGrating,0.1,[25 25 0]); %17 17
else
    GenMaskFromIllu(myGratingTh,0.1,[7 7 0]); %for some of the old version results the 3rd dim was not empty
end

global myillu_mask

if dbg_show
%     cat(4,myillu_mask{:})
    cat(4,myillu_mask{1,:}) %fine grating
    cat(4,myillu_mask{2,:}) %coarse grating
 end


%% Thick-slice blind-SIM
 useCuda=1; enableCuda()
 
% myillu_sumcond={5,10,15};  % necessary for 3D SIM with multi-view deconv. Doesnt work for now
myillu_sumcond=[{5,10,15}; {5,10,15}];
% myillu_sumcond={5,10,15,5,10,15};
% myillu_sumcond={10,20,50};  % necessary for 3D SIM with multi-view deconv. Doesnt work for now
 
myillu=[]; %blind deconv
%psf is rhw 2D in focus PSF and psf3D is the 3D psf
% h=extract(h,[size(h,1) size(h,2) sZ]);
% scaleZ=300; %Try with a bigger scaleZ
% h=kSimPSF({'sX',sX;'sY',sY;'sZ',14;'scaleX',scaleX;'scaleY',scaleY;'scaleZ',scaleZ;'confocal',0;'na',NA;'lambdaEm',lambda;'ri',ri});

% Best results so far
if OldVersion
    DeconvBorders=[0 0 floor(sZ/2)-1]; %Old version
%     [myDeconvBlindT resIllu]=GenericDeconvolution(img,h,[30 5 25 25],'LeastSqr','B',{{'ForcePos',[];'GS',1e-2},{'ForcePos',[]}},[1,1,10],DeconvBorders,[],useCuda,1);
% myillu=[]; %blind deconv
    if withAbberation
        [myDeconvBlindT2 resIllu]=GenericDeconvolution(img,h,[30 5 25 25],'LeastSqr','B',{{'ForcePos',[];'GS',1e-2},{'ForcePos',[]}},[1,1,10],DeconvBorders,[],useCuda,1);
        saveoldD=myDeconvBlindT2;
    else
        [myDeconvBlindT2 resIllu]=GenericDeconvolution(img,h,[30 5 25 25],'LeastSqr','B',{{'ForcePos',[];'GR',9},{'ForcePos',[]}},[1,1,10],DeconvBorders,[],useCuda,1);
        saveold=myDeconvBlindT2;
    end
    
else
    DeconvBorders=[0 0 floor(sZ/2)]; %New version
    if withAbberation
%         [myDeconvBlindT1 resIllu]=GenericDeconvolution(img,h,[30 5 25 5],'LeastSqr','B',{{'ForcePos',[];'GS',0.003},{}},[1,1,20],DeconvBorders,[],useCuda,1);
        
        myillu=[];
        [myDeconvBlindT1 resIllu]=GenericDeconvolution(img,h,[100 5 25 10],'LeastSqr','B',{{'ForcePos',[];'GR',20},{}},[1,1,500],DeconvBorders,[],useCuda,1);
        
    else
    % Quite good (with mask only in first slice):
%     [myDeconvBlindT3 resIllu]=GenericDeconvolution(img,h,[30 5 25 25],'LeastSqr','B',{{'ForcePos',[];'GS',0.005},{}},[1,1,100],DeconvBorders,[],useCuda,1);
    %Better:
%     [myDeconvBlindT3 resIllu]=GenericDeconvolution(img,h,[30 5 25 25],'LeastSqr','B',{{'ForcePos',[];'GR',0.5},{}},[1,1,100],DeconvBorders,[],useCuda,1);
        [myDeconvBlindT3 resIllu]=GenericDeconvolution(img,h,[100 5 25 25],'LeastSqr','B',{{'ForcePos',[];'GR',5},{}},[1,1,100],DeconvBorders,[],useCuda,1);
    %30 and 100 iterations give exactly the same result
    end
end


%% 2D blind-SIM for comparison

if (0)
    
h2D=squeeze(h(:,:,floor(size(h,3)/2))); %PSF2D
DeconvBorders=[0 0];
useCuda=0; disableCuda();
myillu=[]; %blind deconv
if (1) %recalculate the mask for size reason
map=abs(xx([size(a,1) size(a,2) 1])) < 50;
% AbberationMap=-gaussf(map,10)*pi; %That's the one!
AbberationMap=0; 
myGrating2D=GenSIM([size(a,1) size(a,2) 1],SIMDist,3,3,1,AbberationMap); %Illu pattern 2D
AberrationTolerance=10;
myillu_mask=cell(1,numel(myGrating2D));
for v=1:numel(myGrating2D)
    MaskThresh=20;
    myillu_mask{v}=riftshift(dilation(dip_image(abs(rftshift(rft(myGrating2D{v})))> MaskThresh),AberrationTolerance));
end
cat(4,myillu_mask{:})
end

useCuda=1; enableCuda()
[myDeconvBlind2D resIllu]=GenericDeconvolution(img,h2D,[50 5 25 5],'LeastSqr','B',{{'ForcePos',[];'GS',7e-1},{}},[1,1],DeconvBorders,[],useCuda,1);
myillu=[];
% [myDeconvBlind2D1 resIllu]=GenericDeconvolution(img,h2D,[50 5 25 5],'LeastSqr','B',{{'ForcePos',[];'GS',5e-1},{}},[1,1],DeconvBorders,[],useCuda,1);

myillu=[];
% [myDeconvBlind2D2 resIllu]=GenericDeconvolution(img,h2D,[50 5 25 5],'LeastSqr','B',{{'ForcePos',[];'GR',1000},{}},[1,1],DeconvBorders,[],useCuda,1);

end
