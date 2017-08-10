% [img,obj,h,myspeckles,myspecklesIdeal,AbberationMap]=SimImages_Aurelie(sizes,scales,lambda,NA,ri,NumPhotons,NumSpokes,NumSpeck,Zfoc,Zout,Offset,use_obj,use_speckles,use_kSimPSF,sw_AbberationMap)

function [img,obj,h,myspeckles,myspecklesIdeal,AbberationMap]=SimImages_Aurelie(sizes,scales,lambda,NA,ri,NumPhotons,NumSpokes,NumSpeck,Zfoc,Zout,Offset,use_obj,use_speckles,use_kSimPSF,sw_AbberationMap)

%Switches:

if nargin < 15
    sw_AbberationMap=1; % 3 possibilities for the aberration map: (2 implies stronger aberrations on the illumination pattern)
                    %3 = orka figure (only for 2D objects); 0: no
                    %aberration map
end
if nargin < 14
    use_obj='Spoke'; % Choose the object (possibilities are: Spoke, target, pollen, 3dObj, 2dBead, shell, 3dBead)
end
if nargin < 13
    use_kSimPSF=1; %If 1, uses the function kSimPSF to create the psf
end
if nargin < 12
    use_speckles=0; %Choose 1 if the illumination is a speckle pattern and 0 if grating pattern
end
if nargin < 11
    Offset=0;  % 10
end
if nargin < 10
    Zout=6; % Slice index of the out-of-focus information
end
if nargin < 9
    Zfoc=3; % Slice index of the in-focus information
end
if nargin < 8
    NumSpeck=80;  % Number of Speckle images (used only if use_speckles=1)
end
if nargin < 7
    NumSpokes=40;
end
if nargin < 5
    ri=1.52; %Refractive index
end
if nargin < 4
    NA=1.49; %Numerical aperture
end
if nargin < 3
    lambda=500;
end
if nargin < 2
    scales=[lambda/20,lambda/20,200];   % pixel sizes in nm
end
if nargin < 1
    sizes=[200,160,9]; %sZ: number of slices (10 is good for Spokes object, 20 is good for beadSample)
end
scaleX=scales(1);
scaleY=scales(2);
scaleZ=scales(3); 
sX=sizes(1);sY=sizes(2);sZ=sizes(3); %sZ: number of slices (10 is good for Spokes object, 20 is good for beadSample)
if nargin < 6
    NumPhotons=100000 * (4*NA)^2/scaleX/scaleY;    % Photon number in Nyquist pixels
end
    
%% Does not work for thick samples (see speckleTest_TickSlice)
%%
global cuda_enabled;
if ~isempty(cuda_enabled)
    disableCuda();
end

%% User-defined parameters:
%NumPhotons=0;
% SIMDist=lambda/(1.4*NA)/scaleX;  % Emerics Simulation
SIMDist=lambda/(1.8*NA)/scaleX; % at 90% to cutoff
% SIMDist=lambda/(2.8*NA)/scaleX;  % higher Frequency (beyond cutoff)


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
    myspeckles=GenSpeckles(size(obj),0.1,2,NumSpeck); myspecklesIdeal=myspeckles;
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
end

%%
otfs=cell(1,1);
img=cell(numel(myspeckles),1);
otfs{1} = ft(h);
normFac=1;
if NumPhotons > 0   % calcualte normalization only once as speckles may change
    normFac=NumPhotons;
end
obj=Offset+normFac*obj/max(obj);

for p=1:numel(myspeckles)
    fobj = ft(obj.*myspeckles{p});
    mcconv=sqrt(prod(size(obj))) * real(ift(fobj .* otfs{1}));
    mcconv(mcconv<0)=0;
    if NumPhotons > 0
        img{p}=noise(mcconv,'poisson');  % put some noise on the image
    else
        img{p}=mcconv;  % put no noise
    end
end
