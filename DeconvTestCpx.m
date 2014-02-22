% DeconvTestCpx : Test a simple deconvolution using a complex valued PSF rather than a real valued PSF. This is for example used in MRI
NumPhotons=10;
NumPhotons=10000;
Offset=0;  % 10

a=readim('chromo3d');
hw=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',40;'scaleY',40;'scaleZ',100;'confocal',0});
hc=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',40;'scaleY',40;'scaleZ',100;'confocal',1});
hc=hc/max(hc)*max(hw)*0.9;  % roughly realistic
h = hw+i*hc;   % A complex valued PSF
obj=a-i*a/2;   % A complex valued dummy object


fobj = ft(obj);
otfs = ft(h);
mcconv=sqrt(prod(size(obj))) * (ift(fobj .* otfs));  % wii
img=noise(real(mcconv),'gaussian',10)+i*noise(imag(mcconv),'gaussian',10);  % put some noise on the image
%mcconv=norm3d*real(ift(ft(a) .* ft(h)));

%
useCuda=0;
myDeconv=GenericDeconvolution(img,h,15,'LeastSqr',[],{'Complex',[]},[1 1 1],[0 0 0],[],useCuda); st=cat(1,img{1},myDeconv)
    

