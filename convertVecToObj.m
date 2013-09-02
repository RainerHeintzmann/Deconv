function out=convertVecToObj(myinput,DataSize)
PSize=DataSize;PSize(2)=DataSize(1);PSize(1)=DataSize(2);  % To avoid the transpose
global ComplexObj;
if ComplexObj
    out=dip_image(reshape(complex(myinput(1:end/2), myinput(end/2+1:end)),PSize),'scomplex'); % pack two reals to one complex
else
    out=dip_image(reshape(myinput,PSize),'single');  % The reconstruction can be up to 3D, whereas the data might be 4D
end

%if isreal(myinput)  % complex input is not allowed for the minFunc routine (isLegal returns 0)
%else
%    out=dip_image(reshape(myinput,PSize),'scomplex');  % The reconstruction can be up to 3D, whereas the data might be 4D
%end
