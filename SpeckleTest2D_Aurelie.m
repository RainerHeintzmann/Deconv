%% Does not work for thick samples (see speckleTest_TickSlice)
%%
global cuda_enabled;
if ~isempty(cuda_enabled)
    disableCuda();
end

%% User-defined parameters:
%NumPhotons=0;
Offset=0;  % 10
sX=200;sY=160;sZ=1; %sZ: number of slices (10 is good for Spokes object, 20 is good for beadSample)
lambda=500;
scaleX=lambda/20;   % pixel sizes in nm
scaleY=lambda/20;
scaleZ=200; 
NA=1.49; %Numerical aperture
ri=1.52; %Refractive index
NumSpokes=40; %For the Spoke object. 40
% SIMDist=lambda/(1.4*NA)/scaleX;  % Emerics Simulation
SIMDist=lambda/(1.8*NA)/scaleX; % at 90% to cutoff
% SIMDist=lambda/(2.8*NA)/scaleX;  % higher Frequency (beyond cutoff)
NumPhotons=100000 * (4*NA)^2/scaleX/scaleY;    % Photon number in Nyquist pixels
NumSpeck=80;  % Number of Speckle images (used only if use_speckles=1)
AberrationTolerance=5; % defines the radius of the areas of the illumination mask

Zfoc=3; % Slice index of the in-focus information
Zout=6; % Slice index of the out-of-focus information

% Parameters for the blind SIM deconvolution:
Norm='LeastSqr'; SearchMethod='B';
Lambda=0.01; 
Pen='GR'; 
Neg='NegSqr'; 
DeconvBorders=[0 0];
extStack=1; NumIter=25; 
% Pen='NONE'; 
% Neg='NONE';

% For the metadata file:
date='130301';
SimNum=3; %Simulation number (don't forget to update! Or data gets overwritten)

%Switches:
use_obj='Spoke'; % Choose the object (possibilities are: Spoke, target, pollen, 3dObj, 2dBead, shell, 3dBead)
use_kSimPSF=1; %If 1, uses the function kSimPSF to create the psf
use_speckles=0; %Choose 1 if the illumination is a speckle pattern and 0 if grating pattern
useCuda=0; %Choose 1 if the cuda toolbox is installed
use_sumcond=1; %Choose 1  to introduce a sum condition for the blind deconvolution
use_mask=1; %Choose 1 to introduce an illumination mask for the blind deconvolution
do_comp=1; %Compares the blind deconvolution with other modes (if 0, performs only blind deconvolution)
sw_AbberationMap=1; % 3 possibilities for the aberration map: (2 implies stronger aberrations on the illumination pattern)
                    %3 = orka figure (only for 2D objects); 0: no
                    %aberration map
save_res=0; %Choose 1 to save automatically the images, as well as a textfile containing the metadata

% clear myillu_mask; %in case the global variable was already defined
% clear myillu_sumcond; %in case the global variable was already defined
% % Does not work!

%% Create object and PSF for the simulation

switch use_obj
    case 'Spoke'
        a=GenObj('obj','Spoke','sv',[sX sY sZ], 'NumSpokes', NumSpokes,'Zfoc',Zfoc);
        %Other options: Zout (def6), continuous (def1) makeRound (def0)
        
    case 'target'
        a=readim('Y:\MATLAB\images\resolution_fine');
        % The image is approximately 10* bigger in each direction, so the
        % voxel size for the psf is made also 10* bigger (Aurelie 121122)
        scaleX=scaleX/10;scaleY=scaleY/10;scaleZ=scaleZ/10;
        
    case 'pollen'
        a=GenObj('obj','pollen','sv',[sX sY sZ]);
        %Other options: dphi,dtheta (default 0). Default size [128 128 128]
        
    case '3dObj' %3D object made of a hollow sphere, a line...
        a=GenObj('obj','3dObj','sv',[sX sY sZ],'ScaleX',scaleX,'ScaleY',scaleY,'ScaleZ',scaleZ);
        %Other options: line_width (def2),shell_width (def1)
        
    case '2dBead' %1 Slice bead object (randomly generated positions of beads)
        a=GenObj('obj','2dBead','sv',[sX sY]);
        % Other options:nb (def25),mySigma(def2),strength(def255)
        
    case 'shell' %tilted shell
        a=GenObj('obj','shell','sv',[sX sY sZ],'p',0.3,'dtheta',pi/3);
        % other options: dphi(def0), dpsi(def0),shell_width(def1)
        
    case '3dBead' %3D bead sample with defined planes
        a=GenObj('obj','3dBead','sv',[sX sY sZ],'rand',0);
%     Other options: mySigma, strength, rand, nb
    
end
obj=a; %Rename: this is our object

% Simulate PSF
if ~use_kSimPSF
    k0=2*pi/lambda;
    coord=k0*(rr(sX,sY)*scaleX);
    h=(besselj(1,NA*coord)/coord).^2 * k0^2/pi *scaleX*scaleY;
    h(rr(h)==0)= (0.5*NA)^2 * k0^2/pi *scaleX*scaleY;
    h=h/sum(h);
else
    h=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',scaleX;'scaleY',scaleY;'scaleZ',scaleZ;'confocal',0;'na',NA;'lambdaEm',lambda;'ri',ri});
end
%h=readim('psf.ics');
    % h=h/sum(h);


%% Illumination pattern:

if use_speckles
    % NumSpeck=80;  % Number of Speckle images
    myspeckles=GenSpeckles(size(obj),0.1,2,NumSpeck);
else
 %myspeckles=GenSIM(size(a),SIMDist,3,3,1,0);  % coarse grating, phases and directions and contrast (1.0 means perfect)
    if sw_AbberationMap==1
        AbberationMap=- gaussf(abs(xx(size(a))) < 30,8) *pi;
    elseif sw_AbberationMap==2
        AbberationMap= (xx(obj)*xx(obj))/(100*100) *pi;
    elseif sw_AbberationMap==3
        AbberationMap= gaussf(extract(readim('orka.tif')/250,size(obj)),5) * 2* pi; %Only for 2D objects!
    elseif sw_AbberationMap==0
        AbberationMap=0;
    end
    
    myspeckles=GenSIM(size(a),SIMDist,3,3,1,0,AbberationMap);  % coarse grating, phases and directions and contrast (1.0 means perfect)
    myspecklesIdeal=GenSIM(size(a),SIMDist,3,3,1,0);  % coarse grating, phases and directions and contrast (1.0 means perfect)
end %Aurelie 121113: myspecklesIdeal does not exist if use_speckle=1

global myillu;
myillu=myspeckles;
% cat(4,myillu{:}); % to display

%%
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
acq_img=cat(4,img{:}) % Display acquired image

%mcconv=norm3d*real(ift(ft(a) .* ft(h)));

%%  Extract only the in-focus slice from a 3D simulation
if sZ> 1
    for n= 1:length(myspeckles)
        img{n}=squeeze(img{n}(:,:,Zfoc)); %In-focus slice index Zfoc (slice 3; slice 6 is here the out-of-focus info)
        myspeckles{n}=squeeze(myspeckles{n}(:,:,Zfoc));
        if ~use_speckles
            myspecklesIdeal{n}=squeeze(myspecklesIdeal{n}(:,:,Zfoc)); %Kai&Aurelie 061112
        end
        %  The squeeze on the myspeckles is ONLY needed for the code below, if you want to do a two-D (non-thick slice) reconstruction
    end
    % h=squeeze(h(:,:,mid));
    obj=squeeze(obj(:,:,Zfoc)); 
end


%% Deconvolutions to compare

% Inputs for the generic deconvolution function:
% res=GenericDeconvolution(image,psf,NumIter,Method,Update,lambda,Prior,NegPrior,betas,borderSizes,Variances,useCuda)

if do_comp
%     NumIter=55;
    % Compares the blind deconv to other cases
% 1) The deconvolution assumes ideal illumination, although it was
% aberrated
%  if ~use_speckles 
%      myillu=myspecklesIdeal;
%      myDeconvCI=GenericDeconvolution(img,h,NumIter,'LeastSqr',[],0.001,'NONE','NegSqr',[1,1,1],[0 0],[],useCuda); 
%      deconv_ideal=cat(4,myDeconvCI,obj,img{1}) %Display result
%  end

% 2) The deconvolution takes into account the aberrated illumination
% (should perform better than previous one)
myillu=myspeckles;
myDeconvC=GenericDeconvolution(img,h,NumIter,'LeastSqr',[],0.001,'NONE','NegSqr',[1,1,1],[0 0],[],useCuda); 
% deconv_aberrated=cat(4,myDeconvC,obj,img{1}) %Display result

% 3) Wide-field
myillu=[];
img2=squeeze(sum(cat(4,img{:}),[],4)); %Wide-field image
myDeconvWF=GenericDeconvolution(img2,h,NumIter,'LeastSqr',[],0.001,'NONE','NegSqr',[1,1,1],[0 0],[],useCuda); 
% st=cat(4,myDeconvWF,myDeconvC,obj,myspeckles{1},img{1}) %Display result

end %end of comparison

%global TheObject;
%TheObject=cuda(obj);

%% Preparation for the blind deconvolution
% In this deconvolution method, there is possibly a sum condition. The sum of the intensities
% of the first frames should be homogeneous.
% Also, there is the option that the deconvolution modifies the value of only a certain number of
% pixels inside a mask (i.e., inside of the mask, perfect agreement is
% assumed

% myspeckles=myspecklesNS; % Non-squeezed version for the thick slice (does
% not work)

% Preparation of the parameters
% Frames for which the sum condition should be fulfilled
 if use_sumcond
     global myillu_sumcond;
     myillu_sumcond={3,6,9};  % these are the ones for which the sum condition is each fulfilled
%myillu_sumcond={NumSpeck};  % these are the ones for which the sum condition is each fulfilled
 end
myillu=[];

% Illumination mask
if use_mask
    global myillu_mask;   % confines the variable only to a subspace of illu
    myillu_mask=[]; 
    myillu_mask=cell(1,numel(myspeckles));
% AberrationTolerance=5;
%dip_setboundary('periodic')  % Establishes Periodic Boudary conditions.
for v=1:numel(myspeckles)
    MaskThresh=5;
    %if (useCuda)
        %myillu_mask{v}=riftshift(cuda(dilation(dip_image_force(abs(rftshift(rft(cuda(myspecklesIdeal{v}))))> MaskThresh),AberrationTolerance)));
    %myillu_mask{v}=riftshift(dilation(dip_image_force(abs(rftshift(rft(myspecklesIdeal{v})))> MaskThresh),AberrationTolerance));
    myillu_mask{v}=riftshift(dilation(dip_image(abs(rftshift(rft(myspecklesIdeal{v})))> MaskThresh),AberrationTolerance));
    %else
    %    myillu_mask{v}=dilation(abs(ft(myspecklesIdeal{v}))> MaskThresh,AberrationTolerance);
    % end
end
end

%% Blind deconvolution
%Update: B stands for blind
% 9 total slices and keep extended region at the end

%global EvolIllu; %Aurelie 26022013
%global EvolObj; %Aurelie 26022013

[myDeconvBlind,ResIllu,EvolIllu,EvolObj]=GenericDeconvolution(img,h,NumIter,Norm,SearchMethod,Lambda,Pen,Neg,[1,1,1],DeconvBorders,[],useCuda,extStack); 
resIllu=cat(4,myillu{:});
[myDeconvBlind2,ResIllu,EvolIllu,EvolObj]=GenericDeconvolution(img,h,NumIter,Norm,SearchMethod,Lambda,Pen,Neg,[1,1,1],DeconvBorders,[],useCuda,extStack); 
% [myDeconvBlind resIllu]=GenericDeconvolutionAJ(img,h,NumIter,Norm,SearchMethod,Lambda,Pen,Neg,[1,1,1],DeconvBorders,[],useCuda,extStack); 
% st=cat(4,myDeconvWF/max(myDeconvWF),myDeconvBlind/max(myDeconvBlind),myDeconvC/max(myDeconvC),obj/max(obj),myspeckles{1}/max(myspeckles{1}),img{1}/max(img{1}))

% Plots evolution of the cost functional (Aurelie 26022013)
figure(40);
% plot(EvolIllu);
semilogy(EvolIllu);
title('Evolution of the illumination cost functional');
xlabel('Iteration number');
ylabel ('Value (log scale)');
figure(50);
semilogy(EvolObj);
% plot(EvolObj);
title('Evolution of the object cost functional');
xlabel('Iteration number');
ylabel ('Value (log scale)');

% Visualize the results of all the deconvolution methods in one stack
% results_stack=cat(4,myDeconvWF/numel(img),myDeconvCI,myDeconvC,myDeconvBlind,obj,myspeckles{4},dip_image_force(myillu{4}),img{1})
if ~use_speckles
    if useCuda %can use the function dip_image_force
        results_stack=cat(4,myDeconvWF/numel(img),myDeconvCI,myDeconvC,myDeconvBlind,obj,myspeckles{4},dip_image_force(myillu{4}),img{1});
    else
%         results_stack=cat(4,myDeconvWF/numel(img),myDeconvCI,myDeconvC,myDeconvBlind,obj,myspeckles{4},dip_image(myillu{4}),img{1}); %Aurelie 121112: dip_image_force is for Cuda toolbox
    end
    name1='WFdeconv;CIdeconv;Cdeconv;BlindDeconv;obj;speckle;illu;img';
else %myDeconvCI does not exist (mySpecklesIdeal does not exist)
    if useCuda %can use the function dip_image_force
        results_stack=cat(4,myDeconvWF/numel(img),myDeconvC,myDeconvBlind,obj,myspeckles{4},dip_image_force(myillu{4}),img{1});
    else
        results_stack=cat(4,myDeconvWF/numel(img),myDeconvC,myDeconvBlind,obj,myspeckles{4},dip_image(myillu{4}),img{1});
    end
    name1='WFdeconv;CIdeconv;Cdeconv;BlindDeconv;obj;speckle;illu;img';


dipshow(results_stack,'name',name1) %gives the name as a string to the image
end

%% Compare the quality of the deconvolution with the different methods.

if use_obj=='Spoke'
    % Criterium = contrast of the reconstructed object (spokes) with respect to
    % the distance from the center (far away are coarser structures).
    figure(30); % The radial contrast curve will get the handle 30 (in order to be able to save it)
    RadialContrast(obj,NumSpokes/2,'r'); %Original object
    hold on;
    if do_comp
    RadialContrast(myDeconvC,NumSpokes/2,'b'); % Deconvolution with known aberrated illumination pattern
    RadialContrast(myDeconvWF,NumSpokes/2,'g'); % Wide-field
    end
    RadialContrast(myDeconvBlind,NumSpokes/2,'m'); % Blind deconvolution
%     axis([0 100 -.1 1.5]) % for 40 spokes
    axis([0 50 -.1 1.7]) % for 20 spokes
    legend('Object,','Known Illum','WF-Deconv','Blind estimation');
end

%% Save results

if save_res 
    %Create directory in the labbook
    dir=mkdir('C:\Users\Aurelie\Documents\labbook\BlindSIM',date);
    if dir~=1
        fprintf('Impossible to create directory')
    end
    ResPath=strcat('C:\Users\Aurelie\Documents\labbook\BlindSIM\',date);
%     ext=strcat('simu_',int2str(SimNum));
%     dir=mkdir('ResPath',ext);
%     if dir~=1
%         fprintf('Impossible to create directory')
%     end
%     ResPath=strcat(ResPath,'\',ext);
    %metadata file
    filename=sprintf('metadata_simu_%d',SimNum);
    fid=fopen(strcat(ResPath,'\',filename,'.txt'),'a'); %file ID
    if fid==-1
        fprintf('Impossible to create txt file')
    end
    fprintf(fid, '%s', date);
    fprintf(fid, '\nMetadata for the simulation number %d\n\r', SimNum);
    if use_obj=='Spoke'
        fprintf(fid, '\nSimulated object is a spokeobject with %f spokes. \n\r', NumSpokes);
    end
    fprintf(fid, '\nNumber of photons %f\n\r', NumPhotons);    
    fprintf(fid, '\nSize %f x %f\n\r', sX, sY);
    fprintf(fid, '\nnumber of slices %f \n\r', sZ);
    fprintf(fid, '\nOffset %f\n\r', Offset);
    if use_speckles
        fprintf(fid, '\nIllumination pattern %f speckles \n\r', NumSpeck);
    else
        fprintf(fid, '\nIllumination pattern: grating \n\r');
        fprintf(fid, '\nSize of grating %f\n\r', SIMDist*scaleX);
    end
    fprintf(fid,'\nAbberation map number %d \n\r', sw_AbberationMap);
    
    fprintf(fid, '\nDeconvolution:\n\r');
    fprintf(fid, '\n%d iterations on object and then illumination\r\n', NumIter);
%     fprintf(fid, 'Border sizes %f and is keepextendedstack %d\r\n', DeconvBorders, extStack);
    fprintf(fid, '\nNorm (method - function to optimize) is %s\r\n', Norm);
    fprintf(fid, '\nand search method is %s\r\n', SearchMethod);
    fprintf(fid, '\nwith a penalty of type %s and of Lambda %f\r\n', Pen, Lambda);
    fprintf(fid, '\nNegative penalty %s\r\n', Neg);
    if use_mask
        fprintf(fid,'\nWith illumination mask of radius (in A.U) %f\n\r', AberrationTolerance);
    else
        fprintf(fid,'\nWithout illumination mask \n\r');
    end
    if use_sumcond
        fprintf(fid,'\nWith sum condition \r\n');
    else
        fprintf(fid,'\nWithout sum condition \r\n');
    end
       
    fclose(fid);
    
    %save results
    writeim(myDeconvBlind,strcat(ResPath,'\BlindDeconv'));
    writeim(obj,strcat(ResPath,'\Object'));
    if do_comp
    writeim(myDeconvC,strcat(ResPath,'\Deconv'));
    writeim(myDeconvWF/numel(img),strcat(ResPath,'\WFDeconv'));
    end
    % save the plots:
    if use_obj=='Spoke'
        saveas(30,strcat(ResPath,'\RadialContrast')); 
    end
    saveas(40,strcat(ResPath,'\EvolIllu'));
    saveas(50,strcat(ResPath,'\EvolObj'));
end


%%
if 0 % Not run this part for now
    
myillu=[];
wf=squeeze(sum(cat(4,img{:}),[],4));
myDeconvWF=GenericDeconvolution(wf,h,55,'LeastSqr',[],0,'NONE','NONE',[1,1,1],[0 0 0],[],useCuda); 
%Aurelie 121109: previous line: replace the last 1 by "useCuda"
st=cat(4,img{1},myDeconvWF,myDeconvC,obj,myspeckles{1})


myDeconvGP=GenericDeconvolution(img,h,85,'Poisson',[],0.01,'GR','NegSqr',[1,1,1],[0 0 0],[],0); gtp=cat(1,img{1},myDeconvGP)
myDeconvG=GenericDeconvolution(img,h,85,'LeastSqr',[],0.02,'GR','NegSqr',[1,1,1],[0 0 0],[],0); gt=cat(1,img{1},myDeconvG)
myDeconvTVP=GenericDeconvolution(img,h,85,'Poisson',[],0.01,'TV','NegSqr',[1,1,1],[0 0 0],[],0); tvp=cat(1,img{1},myDeconvTVP)
myDeconvTV=GenericDeconvolution(img,h,85,'LeastSqr',[],0.01,'TV','NegSqr',[1,1,1],[0 0 0],[],0); tv=cat(1,img{1},myDeconvTV)
myDeconvArP=GenericDeconvolution(img,h,85,'Poisson',[],0.05,'AR','NegSqr',[1,1,1],[0 0 0],[],0); arp=cat(1,img{1},myDeconvArP)
myDeconvAr=GenericDeconvolution(img,h,85,'LeastSqr',[],0.2,'AR','NegSqr',[1,1,1],[0 0 0],[],0); ar=cat(1,img{1},myDeconvAr)
%myDeconvGR=GenericDeconvolution(img,h,85,'LeastSqr',0.2,'GR','NegSqr',[1,1,1],[0 0 0],[],0); gr=cat(1,img{1},myDeconvGR)
% other update schemes (Ritchardson Lucy)
myDeconvK=GenericDeconvolution(img,h,15,'Poisson','K',0,'NONE','NegSqr',[1,1,1],[0 0 0],[],0); stp=cat(1,imgv,myDeconvK)
myDeconvRL=GenericDeconvolution(img,h,15,'Poisson','RL',0,'NONE','NegSqr',[1,1,1],[0 0 0],[],0); stp=cat(1,img{1},myDeconvRL)
myDeconvRLL=GenericDeconvolution(img,h,15,'Poisson','RLL',0,'NONE','NegSqr',[1,1,1],[0 0 0],[],0); stp=cat(1,img{1},myDeconvRLL)
myDeconvP=GenericDeconvolution(img,h,15,'Poisson',[],1.0,'NONE','NegSqr',[1,1,1],[0 0 0],[],0); stp=cat(1,img{1},myDeconvP)

cat(4,myDeconvRL,myDeconvK,myDeconvRLL,myDeconvP)

mySingleAr=GenericDeconvolution(img(:,:,:,0),h(:,:,:,0),15,'LeastSqr',0.2,'AR','NegSqr',[1,1,1],0,[],0); ars=cat(1,img(:,:,:,0),myDeconvAr,mySingleAr)


% tiffwrite('Results\Chromo_img.tif',img,'no')
% tiffwrite('Results\Chromo_psf.tif',h,'yes')
% tiffwrite('Results\Chromo_DeconvP15i.tif',myDeconvP,'no')
% tiffwrite('Results\Chromo_Deconv15i.tif',myDeconv,'no')
% tiffwrite('Results\Chromo_DeconvTV.tif',myDeconvTV,'no')
% tiffwrite('Results\Chromo_DeconvTVP.tif',myDeconvTVP,'no')
% tiffwrite('Results\Chromo_DeconvAr.tif',myDeconvAr,'no')
% tiffwrite('Results\Chromo_DeconvArP.tif',myDeconvArP,'no')
% tiffwrite('Results\Chromo_DeconvGR.tif',myDeconvG,'no')
% tiffwrite('Results\Chromo_DeconvGRP.tif',myDeconvGRP,'no')
% aslice=10;
% tiffwrite('Results\ChromoSlice_img.tif',img(:,:,aslice),'yes')
% tiffwrite('Results\ChromoSlice_obj.tif',obj(:,:,aslice),'yes')
% tiffwrite('Results\ChromoSlice_DeconvP15i.tif',myDeconvP(:,:,aslice),'yes')
% tiffwrite('Results\ChromoSlice_Deconv15i.tif',myDeconv(:,:,aslice),'yes')
% tiffwrite('Results\ChromoSlice_DeconvGR.tif',myDeconvG(:,:,aslice),'yes')
% tiffwrite('Results\ChromoSlice_DeconvTV.tif',myDeconvTV(:,:,aslice),'yes')
% tiffwrite('Results\ChromoSlice_DeconvTVP.tif',myDeconvTVP(:,:,aslice),'yes')
% tiffwrite('Results\ChromoSlice_DeconvArGamma0.001.tif',myDeconvAr(:,:,aslice),'yes')
% tiffwrite('Results\ChromoSlice_DeconvArPGamma0.001.tif',myDeconvArP(:,:,aslice),'yes')

myDeconvP=GenericDeconvolution(img,h,15,'Poisson',-1,'NONE','NegSqr',[1,1,1],0,[],0); st=cat(1,img,myDeconvP)
myDeconv=GenericDeconvolution(img,h,15,'LeastSqr',-1,'NONE','NegSqr',[1,1,1],0,[],0); st=cat(1,img,myDeconv)
myDeconvTVP=GenericDeconvolution(img,h,50,'Poisson',0.005,'TV','NegSqr',[1,1,1],0,[],0); tv=cat(1,img,myDeconvTVP)
myDeconvTV=GenericDeconvolution(img,h,50,'LeastSqr',0.2,'TV','NegSqr',[1,1,1],0,[],0); tv=cat(1,img,myDeconvTV)
myDeconvAr=GenericDeconvolution(img,h,50,'LeastSqr',30,'AR','NegSqr',[1,1,1],0,[],0); ar=cat(1,img,myDeconvAr)
myDeconvArP=GenericDeconvolution(img,h,50,'Poisson',10,'AR','NegSqr',[1,1,1],0,[],0); arp=cat(1,img,myDeconvArP)
tiffwrite('Results\img.tif',img,'no')
tiffwrite('Results\psf.tif',h,'yes')
tiffwrite('Results\DeconvP15i.tif',myDeconvP,'no')
tiffwrite('Results\Deconv15i.tif',myDeconv,'no')
tiffwrite('Results\DeconvTV.tif',myDeconvTV,'no')
tiffwrite('Results\DeconvTVP.tif',myDeconvTVP,'no')
tiffwrite('Results\DeconvAr.tif',myDeconvAr,'no')
tiffwrite('Results\DeconvArP.tif',myDeconvArP,'no')


showall=cat(4,img,myDeconv,myDeconvP,myDeconvTV,myDeconvTVP,myDeconvAr,myDeconvArP);
cat(1,repmat(obj/max(obj),[1 1 1 size(showall,4)]),showall/max(showall,[],[1 2 3]))

myDeconv=GenericDeconvolution(img,h,15,'Poisson',-1,'NONE','NegSqr',[1,1,1],0,[],0); st=cat(3,img,myDeconv)
myDeconvTV=GenericDeconvolution(img,h,50,'Poisson',0.005,'TV','NegSqr',[1,1,1],0,[],0); tv=cat(3,img,myDeconvTV)
myDeconvAr=GenericDeconvolution(img,h,50,'LeastSqr',10,'AR','NegSqr',[1,1,1],0,[],0); ar=cat(3,img,myDeconvAr)
myDeconvArP=GenericDeconvolution(img,h,50,'Poisson',50,'AR','NegSqr',[1,1,1],0,[],0); arp=cat(3,img,myDeconvArP)

cat(4,img,myDeconv,myDeconvTV,myDeconvAr,myDeconvArP)

end