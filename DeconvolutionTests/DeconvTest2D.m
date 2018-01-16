NumPhotons=40;
% NumPhotons=10000;
Offset=0;  % 10

% a=readim('Y:\MATLAB\images\resolution_fine.tif');
a=readim('C:\Y\MATLAB\images\resolution_coarse.tif');
h=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',40;'scaleY',40;'scaleZ',100;'confocal',0});
obj=a;

fobj = ft(obj);
otf = ft(h);
mcconv=sqrt(prod(size(obj))) * real(ift(fobj .* otf));
img=noise(Offset+NumPhotons*mcconv/max(mcconv),'poisson');  % put some noise on the image

% img(100,100)=1000; % To look for trouble appearing
% img(30,150)=1000;
%%
if (1) % For the 3D sample
    useCuda=1;
    myDeconv=GenericDeconvolution(img,h,85,'Poisson',[],{'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); 
    myDeconvGP=GenericDeconvolution(img,h,85,'Poisson',[],{'GR',0.02;'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); 
    myDeconvKev=GenericDeconvolution(img,h,85,'Poisson',[],{'Kevran',[0.5 0.5];'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); 
    myDeconvTV=GenericDeconvolution(img,h,85,'Poisson',[],{'TV',[0.01 2];'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); 

    tic
    % myDeconvMAPPR = kMAPPR({img,h},{'n',85;'Gamma',0.01});
    toc
    % gtp=cat(3,img,myDeconv,myDeconvGP,myDeconvMAPPR,myDeconvKev,myDeconvTV,obj)
    gtp=cat(3,img,myDeconv,myDeconvGP,myDeconvKev,myDeconvTV,obj)
end    


%%
if (0)
    tiffwrite('Obj.tiff',obj,'yes')
    tiffwrite('img.tiff',img,'yes')
    tiffwrite('DeconvNoReg.tiff',myDeconv,'yes')
    tiffwrite('DeconvGR.tiff',myDeconvGP,'yes')
    tiffwrite('DeconvKevran.tiff',myDeconvKev,'yes')
    tiffwrite('DeconvTV.tiff',myDeconvTV,'yes')
    tiffwrite('PSF.tiff',h,'yes')
end
