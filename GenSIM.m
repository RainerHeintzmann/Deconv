function [myillu,myamp]=GenSIM(asize,D,numphases,numdirs,contrast,startdir, AbberationMap)
if nargin < 7
    AbberationMap=0;
end
if nargin < 6
    startdir=0;
end
if nargin < 5
    contrast=1;
end
if nargin < 4
    numdirs=3;
end
if nargin < 3
    numphases=3;
end
startdir=startdir*pi/180;
f=1/(2*D);
for d=0:numdirs-1
    alpha=pi*d/numdirs+startdir;
    kx=2*pi*f*sin(alpha);
    ky=2*pi*f*cos(alpha);
    for s=1:numphases
        phi0=2*pi/numphases*(s-1)/2;
        myamp{d*numphases+s}=exp(i*(kx*xx(asize)+ky*yy(asize)+phi0+AbberationMap))+exp(-i*(kx*xx(asize)+ky*yy(asize)+phi0));
        myillu{d*numphases+s}=1*(1-contrast)+contrast*real(myamp{d*numphases+s} .* conj(myamp{d*numphases+s}))/2;
    end
end
