%% Create an object that demonstrates the optical sectioning power as well as the resolution of blind-SIM
useCuda=0; disableCuda();
addpath('C:\Users\user\Documents\MATLAB\Aurelie')
a=readim('Y:\MATLAB\images\resolution_512.tif');
% a=readim('Y:\PublicFiles\MATLAB\images\resolution_512.tif');
% a=readim('Y:\MATLAB\images\resolution_fine.tif');
% sV=[256 256]; %Goal size
% convX=1/floor(size(a,1)./sV(1)); %conversion of size (zoom to apply) in X
% convY=1/floor(size(a,2)./sV(2)); %conversion of size (zoom to apply) in Y
% % a2=resample(a,[0.5 0.5],[0 0]); %Resample image (factor 2)
% a2=resample(a,[convX convY],[0 0],'linear'); %Resample image. Ugly

% b=rotation(a,pi);
% c=rotation(a,pi/2);
d=rotation(a,3*pi/2);
obj=newim([size(a,1),size(a,2),8]);

zFoc=3; %focal slice
% dist=2; %distance between 2 objects
obj(:,:,zFoc)=a; %obj(:,:,zFoc+1)=b*1.5;
%obj(:,:,zFoc+2)=c*1.5; 
obj(:,:,zFoc+3)=d*1.5;
% obj(:,:,zFoc+dist)=b*1.5;
a=obj; %rename to avoid confusion later
clear b
clear c
clear d

%% PSF
% Parameters

sX=size(obj,1); sY=size(obj,2); sZ=size(obj,3);
scaleX=25; scaleY=25; scaleZ=200; % pixel sizes in nm (along z it is always controversary)
NA=1.3;
lambda=500;
ri=1.52;

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

% SIMDist=lambda/(1.4*NA)/scaleX;  % Emerics Simulation
SIMDist=lambda/(1.8*NA)/scaleX; % at 90% to cutoff
% SIMDist=lambda/(2.8*NA)/scaleX;  % higher Frequency (beyond cutoff)
NumPhotons=100000 * (4*NA)^2/scaleX/scaleY;    % Photon number in Nyquist pixels
% NumPhotons=100 * (4*NA)^2/scaleX/scaleY;    % Photon number in Nyquist pixels
NumPhotons=NumPhotons*100;
Offset=0;

%% SIM images

% mean(h,[],[1 2]); %should be flat
% real(ift(ft(obj).*ft(h))); %check that in the X-Z view, the flattening is visible
% idwf=real(ift(ft(obj).*ft(h2))); %check that in the X-Z view, the flattening is visible
% idwf(:,:,zFoc);

% AbberationMap= (xx(obj)*xx(obj))/(100*100) *pi;
 map=abs(xx(size(a))) < 50;
  AbberationMap=-gaussf(map,10)*pi; %That's the one!
% AbberationMap=-gaussf(abs(xx(size(a))) < 30,8) *pi; %Not enough aberrated?
% AbberationMap=-gaussf(abs(xx(size(a))) < 50) *pi;
AbberationMap=0;
myGrating=GenSIM(size(a),SIMDist,3,3,1,AbberationMap); %Illu pattern
myGratingTh=GenSIM(size(a),SIMDist,3,3,1,0); %Illu pattern non distorted
global myillu;
myillu=myGrating;

otfs=cell(1,1);
img=cell(numel(myillu),1);
otfs{1} = ft(h);
normFac=1;
if NumPhotons > 0   % calcualte normalization only once as speckles may change
    normFac=NumPhotons;
end
obj=Offset+normFac*obj/max(obj);

for p=1:numel(myillu)
    fobj = ft(obj.*myillu{p});
    mcconv=sqrt(prod(size(obj))) * real(ift(fobj .* otfs{1}));
    mcconv(mcconv<0)=0;
    if NumPhotons > 0
        img{p}=noise(mcconv,'poisson');  % put some noise on the image
    else
        img{p}=mcconv;  % put no noise
    end
end
acq_imgF=cat(4,img{:})
wfimg=squeeze(mean(acq_imgF,[],4)); %Wide-field image

%%  Extract only the in-focus slice from a 3D simulation

if sZ> 1
    for n= 1:length(myGrating)
        img{n}=squeeze(img{n}(:,:,zFoc)); %In-focus slice index Zfoc (slice 3; slice 6 is here the out-of-focus info)
%         myGrating{n}=squeeze(myGrating{n}(:,:,Zfoc));
        %  Extracting the in-focus slice for the grating is ONLY needed if you want to do a two-D (non-thick slice) reconstruction
    end
    % h=squeeze(h(:,:,mid));
end
acq_img=cat(4,img{:})

%% WF Deconvolution
if (0)
useCuda=1; enableCuda();
myillu=[];
img2=squeeze(sum(cat(4,img{:}),[],4)); %Wide-field image
psf2D=squeeze(h(:,:,floor(size(h,3)/2)));
myDeconvWF=GenericDeconvolution(img2,psf2D,100,'LeastSqr',[],{'ForcePos',[];'GS',0.01},[1,1],[0,0],[],useCuda)
myDeconv3DWF=GenericDeconvolution(wfimg,h,100,'LeastSqr',[],{'ForcePos',[];'GS',0.01},[1,1],[0,0],[],useCuda)

end

%% Thick slice blind-SIM - preparation
 useCuda=0; disableCuda();
 
global myillu_sumcond;
myillu_sumcond={3,6,9};  % these are the ones for which the sum condition is each fulfilled

% sZ=20; %Consider much less planes
% myGrating2=GenSIM([size(a,1) size(a,2) sZ],SIMDist,3,3,1,0); %Illu pattern
myGrating2=myGrating;
% myGrating2=myGratingTh;

global myillu_mask;   % confines the variable only to a subspace of illu
myillu_mask=[]; 
myillu_mask=cell(1,numel(myGrating2));
AberrationTolerance=5; %5 for the non-distorted case
%dip_setboundary('periodic')  % Establishes Periodic Boudary conditions.
for v=1:numel(myGrating2)
    MaskThresh=20;
    %if (useCuda)
        %myillu_mask{v}=riftshift(cuda(dilation(dip_image_force(abs(rftshift(rft(cuda(myspecklesIdeal{v}))))> MaskThresh),AberrationTolerance)));
    %myillu_mask{v}=riftshift(dilation(dip_image_force(abs(rftshift(rft(myspecklesIdeal{v})))> MaskThresh),AberrationTolerance));
    myillu_mask{v}=riftshift(dilation(dip_image(abs(rftshift(rft(myGrating2{v})))> MaskThresh),AberrationTolerance));
    %else
    %    myillu_mask{v}=dilation(abs(ft(myspecklesIdeal{v}))> MaskThresh,AberrationTolerance);
    % end
end
cat(4,myillu_mask{:})

%% Thick-slice blind-SIM
 useCuda=1; enableCuda()
 
myillu=[]; %blind deconv
%psf is rhw 2D in focus PSF and psf3D is the 3D psf
DeconvBorders=[0 0 floor(sZ/2)-1];
% h=extract(h,[size(h,1) size(h,2) sZ]);
% scaleZ=300; %Try with a bigger scaleZ
% h=kSimPSF({'sX',sX;'sY',sY;'sZ',14;'scaleX',scaleX;'scaleY',scaleY;'scaleZ',scaleZ;'confocal',0;'na',NA;'lambdaEm',lambda;'ri',ri});


[myDeconvBlindT2 resIllu]=GenericDeconvolution(img,h,[30 5 25 25],'LeastSqr','B',{{'ForcePos',[];'GS',1e-2},{'ForcePos',[]}},[1,1,10],DeconvBorders,[],useCuda,1);
% % [myDeconvBlindT2 resIllu]=GenericDeconvolution(img,h,[30 5 25 5],'LeastSqr','B',{{'ForcePos',[];'GS',1e-2},{'ForcePos',[]}},[1,1,10],DeconvBorders,[],useCuda,1);
% 
% myillu=[]; %blind deconv
% [myDeconvBlindT3 resIllu]=GenericDeconvolution(img,h,[30 5 25 25],'LeastSqr','B',{{'ForcePos',[];'GS',1e-2},{'ForcePos',[]}},[1,1,20],DeconvBorders,[],useCuda,1);

% cat(3,squeeze(myDeconvBlindT(:,:,floor(size(myDeconvBlindT,3)/2))),squeeze(myDeconvBlindT2(:,:,floor(size(myDeconvBlindT,3)/2))),squeeze(myDeconvBlindT(:,:,floor(size(myDeconvBlindT3,3)/2))))
% cat(3,squeeze(myDeconvBlindT(:,:,floor(size(myDeconvBlindT,3)/2))),squeeze(myDeconvBlindT3(:,:,floor(size(myDeconvBlindT3,3)/2))))

% [myDeconvBlindT resIllu]=GenericDeconvolution(img,h,[30 5 10 25],'LeastSqr','B',{{'ForcePos',[];'GS',1e-2},{'ForcePos',[]}},[1,1,5],DeconvBorders,[],useCuda,1);
% myillu=[]; %blind deconv
% [myDeconvBlindT2 resIllu]=GenericDeconvolution(img,h,[30 5 10 25],'LeastSqr','B',{{'ForcePos',[];'GS',1e-2},{'ForcePos',[]}},[1,1,20],DeconvBorders,[],useCuda,1);
% 
% myillu=[]; %blind deconv
% [myDeconvBlindT3 resIllu]=GenericDeconvolution(img,h,[30 5 25 25],'LeastSqr','B',{{'ForcePos',[];'GR',100},{'ForcePos',[]}},[1,1,5],DeconvBorders,[],useCuda,1);


%% 2D blind-SIM for comparison

if (0)
    
h2D=squeeze(h(:,:,floor(size(h,3)/2))); %PSF2D
DeconvBorders=[0 0];
useCuda=0; disableCuda();
myillu=[]; %blind deconv
if (1) %recalculate the mask for size reason
map=abs(xx([size(a,1) size(a,2) 1])) < 50;
AbberationMap=-gaussf(map,10)*pi; %That's the one!
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
[myDeconvBlind2D resIllu]=GenericDeconvolution(img,h2D,[30 5 25 5],'LeastSqr','B',{{'ForcePos',[];'GS',7e-1},{'ForcePos',[]}},[1,1],DeconvBorders,[],useCuda,1);

end
