function res=convertVecToCpxObj(myinput)  % See also the version for non-complex objects
global aRecon;   % will be overwritten
global myim;
DataSize=size(myim{1});
PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose

aRecon=dip_image(reshape(complex(myinput(1:end/2), myinput(end/2+1:end)),PSize),'scomplex'); % pack two reals to one complex

res=aRecon;
