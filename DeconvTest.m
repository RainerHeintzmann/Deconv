NumPhotons=10;
NumPhotons=10000;
Offset=0;  % 10

if (0)
    a=readim('chromo3d');
    h=readim('psf.ics');
    % h=h/sum(h);
else
    a=readim('chromo3d');
    %a=readim;
    %h=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',20;'scaleY',20;'scaleZ',100;'confocal',1});
    %h=readim('psfWF.ics');
    hw=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',40;'scaleY',40;'scaleZ',100;'confocal',0});
    hc=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',40;'scaleY',40;'scaleZ',100;'confocal',1});
    if (0)
        h{1}=hc;
    else
        hc=hc/max(hc)*max(hw)*0.9;  % roughly realistic
        h = {hw,hc};
    end
end
obj=a;


fobj = ft(obj);
otfs=cell(2,1);
img=cell(2,1);
for p=1:numel(h)
    otfs{p} = ft(h{p});
    mcconv=sqrt(prod(size(obj))) * real(ift(fobj .* otfs{p}));
    img{p}=noise(Offset+NumPhotons*mcconv/max(mcconv),'poisson');  % put some noise on the image
end
%mcconv=norm3d*real(ift(ft(a) .* ft(h)));

%%
if (1) % For the 3D sample
    useCuda=0;
    myDeconv=GenericDeconvolution(img,h,15,'LeastSqr',[],{},[1 1 1],[0 0 0],[],useCuda); st=cat(1,img{1},myDeconv)
    myDeconvC=GenericDeconvolution(img,h,15,'LeastSqr',[],{'NegSqr',0.1},[1,1,1],[0 0 0],[],useCuda); st=cat(1,img{1},myDeconvC)
    myDeconvC2=GenericDeconvolution(img,h,15,'LeastSqr',[],{'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); st=cat(1,img{1},myDeconvC2)
    st=cat(1,img{1},myDeconv,myDeconvC,myDeconvC2)
    
    myDeconvGP=GenericDeconvolution(img,h,85,'Poisson',[],{'GR',0.01;'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); gtp=cat(1,img{1},myDeconvGP)
    myDeconvG=GenericDeconvolution(img,h,85,'LeastSqr',[],{'GR',0.02;'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); gt=cat(1,img{1},myDeconvG)
    myDeconvTVP=GenericDeconvolution(img,h,85,'Poisson',[],{'TV',[0.001 0.01];'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); tvp=cat(1,img{1},myDeconvTVP)
    myDeconvTV=GenericDeconvolution(img,h,85,'LeastSqr',[],{'TV',[0.001 0];'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); tv=cat(1,img{1},myDeconvTV)
    
    myDeconvArP=GenericDeconvolution(img,h,85,'Poisson',[],{'AR',0.02;'ForcePos',[]},[1,1,1],[0 0 0],[],0); arp=cat(1,img{1},myDeconvArP)
    myDeconvAr=GenericDeconvolution(img,h,85,'LeastSqr',[],{'AR',0.05;'ForcePos',[]},[1,1,1],[0 0 0],[],0); ar=cat(1,img{1},myDeconvAr)
    
    % other update schemes (Ritchardson Lucy)
    myDeconvK=GenericDeconvolution(img,h,15,'Poisson','K',{'NegSqr',0.1},[1,1,1],[0 0 0],[],useCuda); stp=cat(1,img{1},myDeconvK)
    myDeconvRL=GenericDeconvolution(img,h,15,'Poisson','RL',{'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); stp=cat(1,img{1},myDeconvRL)
    myDeconvRLL=GenericDeconvolution(img,h,150,'Poisson','RLL',{'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); stp=cat(1,img{1},dip_image_force(myDeconvRLL))
    myDeconvP=GenericDeconvolution(img,h,15,'Poisson',[],{'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); stp=cat(1,img{1},myDeconvP)
    
    cat(4,myDeconvK,myDeconvRLL,myDeconvP,myDeconvTV,myDeconvG)
    
    mySingleAr=GenericDeconvolution(img(:,:,:,0),h(:,:,:,0),15,'LeastSqr',0.2,'AR','NegSqr',[1,1,1],0,[],0); ars=cat(1,img(:,:,:,0),myDeconvAr,mySingleAr)
    
    %%
    tiffwrite('Results\Chromo_img.tif',img,'no')
    tiffwrite('Results\Chromo_psf.tif',h,'yes')
    tiffwrite('Results\Chromo_DeconvP15i.tif',myDeconvP,'no')
    tiffwrite('Results\Chromo_Deconv15i.tif',myDeconv,'no')
    tiffwrite('Results\Chromo_DeconvTV.tif',myDeconvTV,'no')
    tiffwrite('Results\Chromo_DeconvTVP.tif',myDeconvTVP,'no')
    tiffwrite('Results\Chromo_DeconvAr.tif',myDeconvAr,'no')
    tiffwrite('Results\Chromo_DeconvArP.tif',myDeconvArP,'no')
    tiffwrite('Results\Chromo_DeconvGR.tif',myDeconvG,'no')
    tiffwrite('Results\Chromo_DeconvGRP.tif',myDeconvGRP,'no')
    aslice=10;
    tiffwrite('Results\ChromoSlice_img.tif',img(:,:,aslice),'yes')
    tiffwrite('Results\ChromoSlice_obj.tif',obj(:,:,aslice),'yes')
    tiffwrite('Results\ChromoSlice_DeconvP15i.tif',myDeconvP(:,:,aslice),'yes')
    tiffwrite('Results\ChromoSlice_Deconv15i.tif',myDeconv(:,:,aslice),'yes')
    tiffwrite('Results\ChromoSlice_DeconvGR.tif',myDeconvG(:,:,aslice),'yes')
    tiffwrite('Results\ChromoSlice_DeconvTV.tif',myDeconvTV(:,:,aslice),'yes')
    tiffwrite('Results\ChromoSlice_DeconvTVP.tif',myDeconvTVP(:,:,aslice),'yes')
    tiffwrite('Results\ChromoSlice_DeconvArGamma0.001.tif',myDeconvAr(:,:,aslice),'yes')
    tiffwrite('Results\ChromoSlice_DeconvArPGamma0.001.tif',myDeconvArP(:,:,aslice),'yes')
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

