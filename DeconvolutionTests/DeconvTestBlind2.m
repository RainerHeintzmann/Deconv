% Simulates a 3D image stack with an aberrated PSF and tries to recover the object and PSF by blind deconvolution

% NumPhotons=10;
%disableCuda();
NumPhotons=1000;
Offset=0;  % 10

if (1)
    a=extract(readim('chromo3d')-9,[160 140 32]);
else
    a=readtimeseries('Y:\MATLAB\Toolboxes\UnserLibrary\invpblib-release\Example\DeconvEx\3D\reference.tif');
    a=extract(a,[128 128 128],[128 128 64]);
end
    % h=readim('psf.ics');
    scaleX=50;scaleY=50;scaleZ=120;
    NA=1.4;
    lambda=520;
    ri=1.52;
    h=aberratedPSF(0,scaleX,lambda,NA,{'ri',ri;'lambdaEm',lambda;'na',NA;'scaleZ',scaleZ;'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);});
    % ha=aberratedPSF(1.8,scaleX,lambda,NA,{'ri',ri;'scaleZ',scaleZ;'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);});
    ha=aberratedPSF(1.0,scaleX,lambda,NA,{'ri',ri;'lambdaEm',lambda;'na',NA;'scaleZ',scaleZ;'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);});
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

tiffwrite('Y:\MATLAB\Toolboxes\UnserLibrary\invpblib-release\Example\DeconvEx\3D\ImgNoAbberation.tif',dip_image(uint16(img{1})),'no');
tiffwrite('Y:\MATLAB\Toolboxes\UnserLibrary\invpblib-release\Example\DeconvEx\3D\ImgAbberated.tif',dip_image(uint16(img{1})),'no');

tiffwrite('Y:\MATLAB\Toolboxes\UnserLibrary\invpblib-release\Example\DeconvEx\3D\PSFNoAbberration.tif',h,'yes','uint16');
tiffwrite('Y:\MATLAB\Toolboxes\UnserLibrary\invpblib-release\Example\DeconvEx\3D\PSFAbberated.tif',ha,'yes','uint16');

obj=Offset+obj*NumPhotons/max(mcconv);
%%
if (1) % For the 3D sample
    useCuda=0;
    maxTestDim=50;
    NormFac=1e-18;
    % Test the gradient:
    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[1 0 0 0 -1e-9],'LeastSqr','Blbfgs',{{'StartImg',obj;'NormFac',NormFac;'MaxTestDim',maxTestDim},{'NormFac',NormFac;'MaxTestDim',maxTestDim},{'NormFac',NormFac;'MaxTestDim',maxTestDim}},[1,1,1],[0 0 0],[],useCuda);

    % Now try the iteration of the PSF with a known object and the aberrated PSF as start value
    [myDeconvRes,dum,resPSF_perf]=GenericDeconvolution(img,ha,[1 0 0 0 10],'LeastSqr','Blbfgs',{{'StartImg',obj},{},{'NormFac',NormFac}},[1,1,1],[0 0 0],[],useCuda); 
    cat(1,ha(:,70,:),resPSF_perf{1}(:,70,:))
    % Now try the iteration of the PSF with a known object and the non-aberrated PSF as start value
    % NOTE: The OTF ist estimated as complex coefficients inside the whole 3D OTF volume
    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[1 0 0 0 150],'LeastSqr','Blbfgs',{{'StartImg',obj},{},{'NormFac',NormFac}},[1,1,1],[0 0 0],[],useCuda); 
    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[1 0 0 0 150],'LeastSqr','Bcg',{{'StartImg',obj},{},{'NormFac',NormFac}},[1,1,1],[0 0 0],[],useCuda); 
    cat(1,ha(:,70,:),resPSF{1}(:,70,:))

    
    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[1 0 0 0 150],'Poisson','BRL',{{'StartImg',obj},{},{'NormFac',1e-9}},[1,1,1],[0 0 0],[],useCuda); 

    [myDeconvRes,dum,resPSF]=GenericDeconvolution(img{1},h,[1 0 0 0 150],'Poisson','Blbfgs',{{'StartImg',obj/sum(obj)},{},{'NormFac',1e-16;'MaxTestDim',maxTestDim}},[1,1,1],[0 0 0],[],useCuda); 

    [resPSF,dum,dum2]=GenericDeconvolution(img{1}/sum(img{1}),obj/sum(obj),[1 150 0 0 0],'Poisson','Blbfgs',{{'StartImg',h;'NormFac',1e-6},{},{'NormFac',1e-6;'MaxTestDim',maxTestDim}},[1,1,1],[0 0 0],[],useCuda); 

    cat(1,h(:,70,:),ha(:,70,:),resPSF{1}(:,70,:))
    cat(1,h(:,70,:),ha(:,70,:),resPSF{1}(:,70,:),resPSF_perf{1}(:,70,:))
    % Exchanging Object and PSF needs a few normalisations:
    [myDeconvRes2,dum,resPSF2]=GenericDeconvolution(img,obj/sum(obj),[1 150 0 0 0],'LeastSqr','Blbfgs',{{'StartImg',h*sum(img{1})},{},{}},[1,1,1],[0 0 0],[],useCuda); 
    cat(1,h(:,70,:),ha(:,70,:),resPSF{1}(:,70,:),myDeconvRes2{1}(:,70,:)/sum(img{1}))

    % Now a proper blind estimation: Unknown object and unknown PSF
    [myDeconvRes3a,dum,resPSF3a]=GenericDeconvolution(img,h,[5 150 0 0 0],'Poisson','RLL',{{},{},{}},[1,1,1],[0 0 0],[],useCuda); 
    [myDeconvRes3,dum,resPSF3]=GenericDeconvolution(img,h,[5 10 50 0 50],'Poisson','BRLL',{{'StartImg',aRecon},{},{}},[1,1,1],[0 0 0],[],useCuda); 
    [myDeconvRes3,dum,resPSF3]=GenericDeconvolution(img,h,[5 0 0 0 50],'LeastSqr','Blbfgs',{{'ForcePos',1;'StartImg',aRecon},{},{}},[1,1,1],[0 0 0],[],useCuda); 
    [myDeconvRes3,dum,resPSF3]=GenericDeconvolution(img,h,[2 150 0 0 50],'LeastSqr','Blbfgs',{{'ForcePos',1;'StartImg',aRecon},{},{}},[1,1,1],[0 0 0],[],useCuda);

    [myDeconvRes3,dum,resPSF3]=GenericDeconvolution(img,h,[15 150 50 0 50],'LeastSqr','Blbfgs',{{'ForcePos',1;},{},{}},[1,1,1],[0 0 0],[],useCuda);
    cat(1,h(:,70,:),ha(:,70,:),resPSF{1}(:,70,:),resPSF3{1}(:,70,:))

%%  This section uses the "ProjPupil" method of PSF estimation. This means only a 2D distrubution of unknown complex values is estimated, which is projected onto a 3D shell, generated by smart interpolation.
    useCuda=0;
    maxTestDim=50;
    NormFac=1e-18;
    % Now estimate only the 2D pupil
    % Test the gradient:
    if (0)
     [myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[1 0 0 0 -1e-3],'LeastSqr','Blbfgs',{{'StartImg',obj;'NormFac',NormFac;'MaxTestDim',maxTestDim},{'NormFac',NormFac;'MaxTestDim',maxTestDim},{'ProjPupil',[lambda,NA,ri],[scaleX scaleY scaleZ];'NormFac',NormFac,[];'MaxTestDim',maxTestDim,[]}},[1,1,1],[0 0 0],[],useCuda);
    end
     % try the deconv
    obj=obj*sum(img{1})/sum(obj);
%     [myDeconvRes,dum,resPSFi]=GenericDeconvolution(img,h,[1 0 0 0 100],'LeastSqr','Blbfgs',{{'StartImg',obj},{},{'ProjPupil',[lambda,NA,ri],[scaleX scaleY scaleZ]}},[1,1,1],[0 0 0],[],useCuda); 
%     cat(1,h(:,70,:),ha(:,70,:),dip_image_force(resPSFi{1}(:,70,:)))

%     [myDeconvRes,dum,resPSF2]=GenericDeconvolution(img,h,[1 5 0 0 50],'LeastSqr','Blbfgs',{{'StartImg',obj},{},{'ProjPupil',[lambda,NA,ri],[scaleX scaleY scaleZ]}},[1,1,1],[0 0 0],[],useCuda); 
%     cat(1,h(:,70,:),ha(:,70,:),resPSF2{1}(:,70,:))
%     
%     [myDeconvRes,dum,resPSF3]=GenericDeconvolution(img,h,[1 5 0 0 50],'LeastSqr','Blbfgs',{{'StartImg',obj;'ForcePos',1},{},{'ProjPupil',[lambda,NA,ri],[scaleX scaleY scaleZ]}},[1,1,1],[0 0 0],[],useCuda); 
%     cat(1,h(:,70,:),ha(:,70,:),resPSF3{1}(:,70,:))
%     
% NA_rec=1.5;
NA_rec=NA;
[myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[6 5 5 0 40],'LeastSqr','Blbfgs',{{'ForcePos',1},{},{'ProjPupil',[lambda,NA_rec,ri],[scaleX scaleY scaleZ]}},[scaleX scaleY scaleZ],[0 0 9],[],useCuda);
%[myDeconvRes,dum,resPSF]=GenericDeconvolution(img,h,[5 5 5 0 40],'LeastSqr','Blbfgs',{{'ForcePos',1},{},{'ProjPupil',[lambda,NA,ri],[scaleX scaleY scaleZ]}},[1,1,1],[0 0 0],[],useCuda); 
    % stp=cat(1,obj,img{1},myDeconvRes)
cat(1,h(:,128,:)/max(h),ha(:,70,:)/max(ha),(resPSF{1}(:,70,:)/max(resPSF{1})))
    cat(1,h(:,70,:)/max(h),ha(:,70,:)/max(ha),dip_image_force(resPSF{1}(:,70,:)/max(resPSF{1})))
    global savedATF
    phase(savedATF)

 %%
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

