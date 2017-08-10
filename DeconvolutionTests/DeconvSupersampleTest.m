NumPhotons=10;
NumPhotons=10000;
scaleX=40;
scaleZ=100;
Offset=0;  % 10

if (1)
    a=readim('chromo3d');
    h=extract(readim('psf.ics'),size(a));
    % h=h/sum(h);
else
    a=readim('chromo3d');
    %a=readim;
    %h=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',20;'scaleY',20;'scaleZ',100;'confocal',1});
    %h=readim('psfWF.ics');
    hw=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',scaleX;'scaleY',scaleX;'scaleZ',scaleZ;'confocal',0});
    hc=kSimPSF({'sX',size(a,1);'sY',size(a,2);'sZ',size(a,3);'scaleX',scaleX;'scaleY',scaleX;'scaleZ',scaleZ;'confocal',1});
    if (1)
        h{1}=hw;
    else
        hc=hc/max(hc)*max(hw)*0.9;  % roughly realistic
        h = {hw,hc};
    end
end
obj=a;

fobj = ft(obj);
otfs=cell(2,1);
img=cell(1,1);
for p=1:numel(h)
    otfs{p} = ft(h{p});
    mcconv=sqrt(prod(size(obj))) * real(ift(fobj .* otfs{p}));
    img{p}=noise(Offset+NumPhotons*mcconv/max(mcconv),'poisson');  % put some noise on the image
end
%mcconv=norm3d*real(ift(ft(a) .* ft(h)));

%%
img{1}=rift(rft_resize(rft(img{1}),0.5));
h=kSimPSF({'sX',size(a,1)/2;'sY',size(a,2)/2;'sZ',size(a,3)/2;'scaleX',scaleX/2;'scaleY',scaleX/2;'scaleZ',scaleZ/2;'confocal',0});

%%
useCuda=0;
myDeconvS=GenericDeconvolution(img,h,195,'Poisson','RLL',{'Resample',2},[1 1 1],[0 0 0],[],useCuda)
myDeconv=GenericDeconvolution(img,h,195,'Poisson','RLL',{},[1 1 1],[0 0 0],[],useCuda); 
st=cat(1,img{1},myDeconv)
st=cat(1,obj*1000,myDeconvS)
