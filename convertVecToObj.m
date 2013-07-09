function out=convertVecToObj(myinput,DataSize)
PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose
out=dip_image(reshape(myinput,PSize),'single');  % The reconstruction can be up to 3D, whereas the data might be 4D
