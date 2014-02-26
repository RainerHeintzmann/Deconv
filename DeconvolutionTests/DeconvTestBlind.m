% Simulates a 3D image stack with an aberrated PSF and tries to recover the object and PSF by blind deconvolution

% NumPhotons=10;
NumPhotons=10000;
Offset=0;  % 10

    a=extract(readim('chromo3d')-9,[160 140 32]);
    % h=readim('psf.ics');
    scaleX=50;scaleY=50;scaleZ=150;
    NA=1.4;
    lambda=520;
    ri=1.52;
    h=aberratedPSF(0,scaleX,lambda,NA,{'ri',ri;'scaleZ',scaleZ;'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);});
    ha=aberratedPSF(1.8,scaleX,lambda,NA,{'ri',ri;'scaleZ',scaleZ;'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);});
    obj=a;

    cat(1,h(:,70,:),ha(:,70,:))

fobj = ft(obj);
otfs=cell(numel(h),1);
img=cell(numel(h),1);
for p=1:numel(h)
    otfs{p} = ft(ha{p});  % use aberrated PSF for simulations
    mcconv=sqrt(prod(size(obj))) * real(ift(fobj .* otfs{p}));
    img{p}=noise(Offset+NumPhotons*mcconv/max(mcconv),'poisson');  % put some noise on the image
end
obj=Offset+obj*NumPhotons/max(mcconv);
%%
if (1) % For the 3D sample
    useCuda=0;
    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[1 0 0 0 150],'LeastSqr','Blbfgs',{{'StartImg',obj},{},{}},[1,1,1],[0 0 0],[],useCuda); 
    cat(1,h(:,70,:),ha(:,70,:),resPSF{1}(:,70,:))
    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[1 0 0 0 150],'LeastSqr','Blbfgs',{{'StartImg',obj},[],{'ProjPupil',[lambda,NA,ri],[scaleX scaleY scaleZ]}},[1,1,1],[0 0 0],[],useCuda); 
    cat(1,h(:,70,:),ha(:,70,:),resPSF{1}(:,70,:))
    
    
    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,ha,[1 0 0 0 150],'LeastSqr','Blbfgs',{{'StartImg',obj},{},{}},[1,1,1],[0 0 0],[],useCuda); 

    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[1 0 0 0 150],'LeastSqr','Blbfgs',{{'StartImg',obj},{},{'ProjPupil',[512,0.8],[100 100 200]}},[1,1,1],[0 0 0],[],useCuda); stp=cat(1,obj,img{1},myDeconvRes)

    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[20 5 5 0 15],'LeastSqr','Blbfgs',{{},{},{}},[1,1,1],[0 0 0],[],useCuda); stp=cat(1,obj,img{1},myDeconvRes)

    myDeconvRLL=GenericDeconvolution(img,h,[50 15 5 0 5],'Poisson','RLL',{{},{},{}},[1,1,1],[0 0 0],[],useCuda); stp=cat(1,obj,img{1},myDeconvRLL)
    myDeconvBRLL=GenericDeconvolution(img,h,[50 15 5 0 5],'Poisson','BRLL',{{},{},{}},[1,1,1],[0 0 0],[],useCuda); stp=cat(1,obj,img{1},myDeconvBRLL)
    
    stp=cat(1,obj,img{1},myDeconvRLL,myDeconvBRLL)
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

