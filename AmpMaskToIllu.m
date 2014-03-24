% realimage=AmpMaskToIllu(ampmask,myvector)  : fills the vector into the mask region and Fourier transforms back to real space

function realimage=AmpMaskToIllu(ampmask,myvector,maskindex)
tmpim=newim(size(ampmask),'scomplex');
tmpim(ampmask)=myvector;   % writes the complex vector into the Fourier-image
if nargin > 2
    global mytmp;
    if iscell(mytmp)
        mytmp{maskindex}=tmpim+0;
    end
end
if (1)
    tmp2=fixRealRFTValues(tmpim)+0;  % Just to guarantee a good result for the inverse RFT, as certain values need to be real
else
    tmp2=tmpim;  %no fixing in blind iterations!
end
realimage=rift(tmp2)+0;         % go to real space (cuda version)
% realimage=ift(tmpim);         % go to real space (cuda version). Aurelie 28.01.2014
