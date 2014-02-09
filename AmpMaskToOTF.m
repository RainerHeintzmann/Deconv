% realimage=AmpMaskToIllu(ampmask,myvector)  : fills the vector into the mask region and Fourier transforms back to real space

function realimage=AmpMaskToOTF(ampmask,myvector)
realimage=newim(ampmask,'scomplex');
realimage(ampmask)=myvector;   % writes the complex vector into the Fourier-image

%tmpim=rift(tmpim);         % go to real space (cuda version)
%realimage=tmpim;
