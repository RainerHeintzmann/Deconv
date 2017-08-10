global cuda_enabled;
if ~isempty(cuda_enabled)
    disableCuda();
end
%NumPhotons=0;
Offset=0;  % 10
sX=200;sY=200;sZ=7;
lambda=500;
scaleX=lambda/20;   % pixel sizes in nm
scaleY=lambda/20;
scaleZ=200; 
NA=1.49;
ri=1.52;
NumSpokes=40;
SIMDist=lambda/(1.4*NA)/scaleX;  % Emerics Simulation
% SIMDist=lambda/(2.8*NA)/scaleX;  % higher Frequency (beyond cutoff)
NumPhotons=1000000 * (4*NA)^2/scaleX/scaleY;    % Photon number in Nyquist pixels

if (1)
    %a=readim('chromo3d');
    %a=readim('lena');
    %a=readim;
    if sZ<= 1
        a=SpokesObject([sX sY],NumSpokes,1,0);
    else
        a=newim(sX,sY,sZ);
        a1=SpokesObject([sX sY],NumSpokes,1,0);
        a2=abs(xx(sX,sY)-30)< 30;
        a3=abs(xx(sX,sY)+30)< 30;
        a(:,:,3)=a1;
        if (0)
            a(:,:,3)=a1+a3*2;
            a(:,:,3)=a1+a3*2;
            a(:,:,6)=a2*6;
        end
    end
    
    if (0)
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
else
    a=readim('chromo3d');
    %a=readim;
    %h=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',20;'scaleY',20;'scaleZ',100;'confocal',1});
    %h=readim('psfWF.ics');
    hw=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',20;'scaleY',20;'scaleZ',100;'confocal',0});
    hc=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',20;'scaleY',20;'scaleZ',100;'confocal',1});
    if (0)
        h{1}=hc;
    else
        hc=hc/max(hc)*max(hw)*0.9;  % roughly realistic
        h = {hw,hc};
    end
end
obj=a;

if(1)
NumSpeck=80;  % Number of Speckle images
myspeckles=GenSpeckles(size(obj),0.1,2,NumSpeck);
myspecklesIdeal=myspeckles;
else
 %myspeckles=GenSIM(size(a),SIMDist,3,3,1,0);  % coarse grating, phases and directions and contrast (1.0 means perfect)
AbberationMap=- gaussf(abs(xx(size(a))) < 30,8) *pi;
% AbberationMap= (xx(obj)*xx(obj))/(100*100) *pi;
% AbberationMap= gaussf(extract(readim('orka.tif')/250,size(obj)),5) * 2* pi;
 
% A two-beam simulation is assumed. Thus the out-of-focus regions get
% modulated just as well
 myspeckles=GenSIM(size(a),SIMDist,3,3,1,0,AbberationMap);  % coarse grating, phases and directions and contrast (1.0 means perfect)
 myspecklesIdeal=GenSIM(size(a),SIMDist,3,3,1,0);  % coarse grating, phases and directions and contrast (1.0 means perfect)
end

global myillu;
myillu=myspeckles;

otfs=cell(1,1);
img=cell(numel(myillu),1);
otfs{1} = ft(h);
normFac=1;
if NumPhotons > 0   % calcualte normaization only once as speckles may change
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
cat(4,img{:})

%mcconv=norm3d*real(ift(ft(a) .* ft(h)));

%%  Extract only the in-focus slice from a 3D simulation, this also needs to be run for the thick-slice recon
for n= 1:length(myspeckles)
    img{n}=squeeze(img{n}(:,:,3));
end
mid=floor(size(h,3)/2);
h=squeeze(h(:,:,mid));
obj=squeeze(obj(:,:,3));
%%  This is ONLY needed for the code below, if you want to do a two-D (non-thick slice) reconstruction
for n= 1:length(myspeckles)
    myspeckles{n}=squeeze(myspeckles{n}(:,:,3));
    myspecklesIdeal{n}=squeeze(myspecklesIdeal{n}(:,:,3));
end

%%
if (0)
    useCuda=1;enableCuda();
else
    useCuda=0; % disableCuda();
end

myillu=[];  % myspecklesIdeal; %Blind
NumIter=[50 5 25 5];

myDeconvCI=GenericDeconvolution(img,h,NumIter,'LeastSqr','B',{{'NegSqr',0.001},{'NegSqr',0.001}},[1,1,1],[0,0],[],useCuda,0); 
% myDeconvCI=GenericDeconvolution(img,h,55,'LeastSqr',[],0.001,'NONE','NegSqr',[1,1,1],[0 0],[],useCuda); 
st=cat(4,myDeconvCI,obj,img{1})

%%
myillu=myspeckles;

myDeconvC=GenericDeconvolution(img,h,55,'LeastSqr',[],0.001,'NONE','NegSqr',[1,1,1],[0 0],[],useCuda); 

%myDeconvC=GenericDeconvolution(img,h,155,'LeastSqr',[],0.01,'GR','NegSqr',[1,1,1],[0 0],[],useCuda); 
%myDeconvC=GenericDeconvolution(img,h,155,'Poisson','RLL',0.0,'NONE','NegSqr',[1,1,1],[0 0],[],useCuda); 
st=cat(4,myDeconvC,myDeconvCI,obj,img{1})

myillu=[];
img2=squeeze(sum(cat(4,img{:}),[],4));
myDeconvWF=GenericDeconvolution(img2,h,55,'LeastSqr',[],0.001,'NONE','NegSqr',[1,1,1],[0 0],[],useCuda); 
% st=cat(4,myDeconvWF,myDeconvC,obj,myspeckles{1},img{1})

%global TheObject;
%TheObject=cuda(obj);

%%
useCuda=1;
global myillu_sumcond;
%myillu_sumcond={3,6,9};  % these are the ones for wich the sum condition is each fulfilled
myillu_sumcond={NumSpeck};  % these are the ones for wich the sum condition is each fulfilled
myillu=[];

%%
global myillu_mask;   % confines the variable only to a subspace of illu
myillu_mask=[];
myillu_mask=cell(1,numel(myspeckles));
AberrationTolerance=5;
%dip_setboundary('periodic')  % Establishes Periodic Boudary conditions.
for v=1:numel(myspeckles)
    MaskThresh=5;
    if (useCuda)
        myillu_mask{v}=riftshift(cuda(dilation(dip_image_force(abs(rftshift(rft(cuda(myspecklesIdeal{v}))))> MaskThresh),AberrationTolerance)));
    else
        myillu_mask{v}=dilation(abs(ft(myspecklesIdeal{v}))> MaskThresh,AberrationTolerance);
    end
end

%%
myillu=[];
% myillu_mask=[];  % here you can chose whether to use the Fourier mask or not
if (0)
    useCuda=1;enableCuda();
else
    useCuda=0;disableCuda();
end
myDeconvBlind=GenericDeconvolution(img,h,25,'LeastSqr','B',0.02,'GR','NegSqr',[1,1,1],[0 0 0],[],useCuda,1); % [0 0 4] for thick slice
%myDeconvBlind=GenericDeconvolution(img,h,25,'Poisson','BRLL',1e-5,'NONE','NegSqr',[1,1,1],[0 0 0],[],useCuda,1); % [0 0 4] for thick slice
%myDeconvBlind=GenericDeconvolution(img,h,55,'LeastSqr','BRLL',0.01,'TV','NegSqr',[1,1,1],[0 0 0],[],useCuda,1); % [0 0 4] for thick slice
% st=cat(4,myDeconvWF/max(myDeconvWF),myDeconvBlind/max(myDeconvBlind),myDeconvC/max(myDeconvC),obj/max(obj),myspeckles{1}/max(myspeckles{1}),img{1}/max(img{1}))

st=cat(4,myDeconvWF/numel(img),myDeconvCI,myDeconvC,myDeconvBlind,obj,myspeckles{4},dip_image_force(myillu{4}),img{1})

%%
figure;
RadialContrast(obj,NumSpokes/2,'r');
hold on;
RadialContrast(myDeconvCI,NumSpokes/2,'m');
RadialContrast(myDeconvC,NumSpokes/2,'b');
RadialContrast(myDeconvWF,NumSpokes/2,'g');
RadialContrast(myDeconvBlind,NumSpokes/2,'c');
axis([0 100 -.1 1.5])
legend('Object,','Known Illumm ideal assumed','Known Illum','WF-Deconv','Blind estimation');

%%
if 1
myillu=[];
wf=squeeze(sum(cat(4,img{:}),[],4));
myDeconvWF=GenericDeconvolution(wf,h,55,'LeastSqr',[],0,'NONE','NONE',[1,1,1],[0 0 0],[],1); 
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
else
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
end

showall=cat(4,img,myDeconv,myDeconvP,myDeconvTV,myDeconvTVP,myDeconvAr,myDeconvArP);
cat(1,repmat(obj/max(obj),[1 1 1 size(showall,4)]),showall/max(showall,[],[1 2 3]))

myDeconv=GenericDeconvolution(img,h,15,'Poisson',-1,'NONE','NegSqr',[1,1,1],0,[],0); st=cat(3,img,myDeconv)
myDeconvTV=GenericDeconvolution(img,h,50,'Poisson',0.005,'TV','NegSqr',[1,1,1],0,[],0); tv=cat(3,img,myDeconvTV)
myDeconvAr=GenericDeconvolution(img,h,50,'LeastSqr',10,'AR','NegSqr',[1,1,1],0,[],0); ar=cat(3,img,myDeconvAr)
myDeconvArP=GenericDeconvolution(img,h,50,'Poisson',50,'AR','NegSqr',[1,1,1],0,[],0); arp=cat(3,img,myDeconvArP)

cat(4,img,myDeconv,myDeconvTV,myDeconvAr,myDeconvArP)

