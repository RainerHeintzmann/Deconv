function [acurve,bincenters]=radialmean2(img,pixelsizes,binsize,minr)
if nargin < 2
    pixelsizes=size(img)*0+1;
end
if nargin < 3
    binsize=1;
end
if nargin < 4
    minr=1;
end
myr2=0;
for n=1:ndims(img)
    myr2=myr2+(pixelsizes(n)*ramp(size(img),n))^2;
end
myr=sqrt(myr2);

nbins=max(myr)/binsize-minr/binsize;

% myr_bin = round((myr-minr)/binsize);

mm=ceil(max(myr));
bincenters=[minr:binsize:mm];
acurve=newim([size(bincenters,2),1],datatype(img));
for n=1:size(bincenters,2)
    mymask=(myr>=bincenters(n)-binsize/2) & (myr<bincenters(n)+binsize/2);
    acurve(n-1)=mean(img,mymask);
end
