% realimage=AmpMaskToIllu(ampmask,myvector)  : fills the vector into the mask region and Fourier transforms back to real space

function realimage=AmpMaskToIllu(ampmask,myvector)
global mytmp;
tmpim=newim(ampmask,'scomplex');
tmpim(ampmask)=myvector;   % writes the complex vector into the Fourier-image
mytmp=tmpim;
tmpim=rift(tmpim);         % go to real space (cuda version)
% tmpim=ift(tmpim);         % go to real space (cuda version). Aurelie 28.01.2014

realimage=tmpim;
