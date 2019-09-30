% realimage=AmpMaskToIllu(ampmask,myvector)  : fills the vector into the mask region and Fourier transforms back to real space

function otf=ReplaceAmpMask(otf,ampmask,myvector)
otf(ampmask)=myvector;   % writes the complex vector into the Fourier-image

%tmpim=rift(tmpim);         % go to real space (cuda version)
%realimage=tmpim;
