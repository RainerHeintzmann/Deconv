NumPhotons=100;
Offset=10;

a=readim('chromo3d');
origSize=size(a);
a=extract(a-10,origSize*2);
h=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',20;'scaleY',20;'scaleZ',100;'confocal',0});

obj=a;
obj(78,130,15)=max(obj)*100;

mcconv=sqrt(prod(size(a)))*real(ift(ft(obj) .* ft(h)));
mcconv=extract(mcconv,origSize);
hc=extract(h,origSize);

img=noise(Offset+NumPhotons*mcconv/max(mcconv),'poisson');  % put some extra noise on the image

myDeconvTVP=GenericDeconvolution(img,hc,50,'Poisson',0.01,'TV','NegSqr',1,[1,1,1],[3 3 3],0); tv=cat(1,img,myDeconvTVP)

%mcconv=norm3d*real(ift(ft(a) .* ft(h)));

img=noise(Offset+NumPhotons*mcconv/max(mcconv),'poisson');  % put some extra noise on the image

myDeconvP=GenericDeconvolution(img,h,15,'Poisson',0,'NONE','NegSqr',1,[1,1,1],0); st=cat(1,img,myDeconvP)
myDeconv=GenericDeconvolution(img,h,15,'LeastSqr',0,'NONE','NegSqr',1.0,[1,1,1],0); st=cat(1,img,myDeconv)
myDeconvTVP=GenericDeconvolution(img,h,50,'Poisson',0.01,'TV','NegSqr',1,[1,1,1],0); tv=cat(1,img,myDeconvTVP)
myDeconvTV=GenericDeconvolution(img,h,50,'LeastSqr',0.2,'TV','NegSqr',1.0,[1,1,1],0); tv=cat(1,img,myDeconvTV)
myDeconvAr=GenericDeconvolution(img,h,50,'LeastSqr',2,'Arigovindan','NegSqr',1.0,[1,1,1],0); ar=cat(1,img,myDeconvAr)
myDeconvArP=GenericDeconvolution(img,h,50,'Poisson',0.05,'Arigovindan','NegSqr',1,[1,1,1],0); arp=cat(1,img,myDeconvArP)

tiffwrite('Results\Chromo_img.tif',img,'no')
tiffwrite('Results\Chromo_psf.tif',h,'yes')
tiffwrite('Results\Chromo_DeconvP15i.tif',myDeconvP,'no')
tiffwrite('Results\Chromo_Deconv15i.tif',myDeconv,'no')
tiffwrite('Results\Chromo_DeconvTV.tif',myDeconvTV,'no')
tiffwrite('Results\Chromo_DeconvTVP.tif',myDeconvTVP,'no')
tiffwrite('Results\Chromo_DeconvAr.tif',myDeconvAr,'no')
tiffwrite('Results\Chromo_DeconvArP.tif',myDeconvArP,'no')
slice=10;
tiffwrite('Results\ChromoSlice_img.tif',img(:,:,slice),'yes')
tiffwrite('Results\ChromoSlice_obj.tif',obj(:,:,slice),'yes')
tiffwrite('Results\ChromoSlice_DeconvP15i.tif',myDeconvP(:,:,slice),'yes')
tiffwrite('Results\ChromoSlice_Deconv15i.tif',myDeconv(:,:,slice),'yes')
tiffwrite('Results\ChromoSlice_DeconvTV.tif',myDeconvTV(:,:,slice),'yes')
tiffwrite('Results\ChromoSlice_DeconvTVP.tif',myDeconvTVP(:,:,slice),'yes')
tiffwrite('Results\ChromoSlice_DeconvArGamma0.001.tif',myDeconvAr(:,:,slice),'yes')
tiffwrite('Results\ChromoSlice_DeconvArPGamma0.001.tif',myDeconvArP(:,:,slice),'yes')

cat(4,img,myDeconv,myDeconvP,myDeconvTV,myDeconvTVP,myDeconvAr,myDeconvArP)


myDeconv=GenericDeconvolution(img,h,15,'Poisson',-1,'NONE','NegSqr',1.0,[1,1,1],0); st=cat(3,img,myDeconv)
myDeconvTV=GenericDeconvolution(img,h,50,'Poisson',0.005,'TV','NegSqr',1.0,[1,1,1],0); tv=cat(3,img,myDeconvTV)
myDeconvAr=GenericDeconvolution(img,h,50,'LeastSqr',10,'Arigovindan','NegSqr',1.0,[1,1,1],0); ar=cat(3,img,myDeconvAr)
myDeconvArP=GenericDeconvolution(img,h,50,'Poisson',50,'Arigovindan','NegSqr',1.0,[1,1,1],0); arp=cat(3,img,myDeconvArP)

cat(4,img,myDeconv,myDeconvTV,myDeconvAr,myDeconvArP)

