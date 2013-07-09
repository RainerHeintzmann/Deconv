function res=psf(myim)

res=kSimPSF({'sZ',size(myim,3);'sX',size(myim,1);'sY',size(myim,2)});
