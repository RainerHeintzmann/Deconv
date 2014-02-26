Offset=0;  % 10
sX=200;sY=160;sZ=1; %sZ: number of slices (10 is good for Spokes object, 20 is good for beadSample)
lambda=500;
scaleX=lambda/20;   % pixel sizes in nm
scaleY=lambda/20;
scaleZ=200;
NA=1.49; %Numerical aperture
ri=1.52;
SIMDist=lambda/(1.8*NA)/scaleX; % at 90% to cutoff
NumPhotons=100000 * (4*NA)^2/scaleX/scaleY;    % Photon number in Nyquist pixels
NumSpeck=80;  % Number of Speckle images (used only if use_speckles=1)
AberrationTolerance=5; % defines the radius of the areas of the illumination mask

obj = imread('einstein2.tif');
a = obj;

myillu=GenSIM(size(a),SIMDist,3,3,1,0);

h=kSimPSF({'sX',size(a,1);'sY',size(a,2);'scaleX',scaleX;'scaleY',scaleY;'scaleZ',scaleZ;'confocal',0;'na',NA;'lambdaEm',lambda;'ri',ri});

otfs=cell(1,1);
img=cell(numel(myillu),1);
otfs{1} = ft(h);
normFac=1;

obj=Offset+normFac*(obj./max(obj, 1));
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
acq_img=cat(4,img{:})
