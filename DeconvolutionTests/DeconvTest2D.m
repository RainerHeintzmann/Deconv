NumPhotons=10;
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
    myDeconvGP=GenericDeconvolution(img,h,85,'Poisson',[],{'GR',0.01;'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); 

    myDeconv=GenericDeconvolution(img,h,85,'Poisson',[],{'ForcePos',[]},[1,1,1],[0 0 0],[],useCuda); 

    tic
    myDeconvMAPPR = kMAPPR({img,h},{'n',85;'Gamma',0.01});
    toc
    gtp=cat(3,img,myDeconvGP,myDeconvMAPPR,myDeconv,obj)
end    
    

