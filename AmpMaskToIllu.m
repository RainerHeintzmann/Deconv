% realimage=AmpMaskToIllu(ampmask,myvector)  : fills the vector into the mask region and Fourier transforms back to real space

function realimage=AmpMaskToIllu(ampmask,myvector,maskindex)
tmpim=newim(ampmask,'scomplex');
tmpim(ampmask)=myvector;   % writes the complex vector into the Fourier-image
if nargin > 2
    global mytmp;
    if iscell(mytmp)
        mytmp{maskindex}=tmpim;
    end
end
tmpim=fixRealRFTValues(tmpim);  % Just to guarantee a good result for the inverse RFT, as certain values need to be real
tmpim=rift(tmpim);         % go to real space (cuda version)
% tmpim=ift(tmpim);         % go to real space (cuda version). Aurelie 28.01.2014

realimage=tmpim;
