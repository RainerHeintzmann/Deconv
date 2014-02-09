function res=convertRVecToCpxObj(myinput)  % See also the version for non-complex objects
global aRecon;  % Will be overwritten
global myim;
DataSize=size(myim{1});
PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose

aRecon=dip_image(reshape(complex(myinput),PSize),'scomplex'); % pack two reals to one complex

res=aRecon;