% correctedRFT=fixRealRFTValues(anRFT) : Forces the RFT Values at symmetric positions to be real. This is needed to get correct rft results
function outRFT=fixRealRFTValues(anRFT)
anRFT=anRFT+0;
if ndims(anRFT) == 1
    anRFT(0)=real(anRFT(0));
    mx=floor(size(anRFT,1)/2);
    tmp=real(anRFT(mx));         % All those tmp assignments are needed for Cuda not getting a hickup.
    anRFT(mx)=tmp;
elseif ndims(anRFT) == 2
    tmp=real(anRFT(0,0));
    anRFT(0,0)=tmp;
    tmp=real(anRFT(0,end));
    anRFT(0,end)=tmp;
    mx=floor(size(anRFT,1)/2);
    if 2*mx==size(anRFT,1)
        tmp=real(anRFT(mx,0));
        anRFT(mx,0)=tmp;
        tmp=real(anRFT(mx,end));
        anRFT(mx,end)=tmp;
    end
elseif ndims(anRFT) == 3
    mx=floor(size(anRFT,1)/2);
    tmp=real(anRFT(0,0,0));
    anRFT(0,0,0)=tmp;
    tmp=real(anRFT(0,end,0));
    anRFT(0,end,0)=tmp;
    if 2*mx==size(anRFT,1)
        tmp=real(anRFT(mx,0,0));
        anRFT(mx,0,0)=tmp;
        tmp=real(anRFT(mx,end,0));
        anRFT(mx,end,0)=tmp;
    end
    mz=floor(size(anRFT,3)/2);
    if 2*mz==size(anRFT,3)
        tmp=real(anRFT(0,0,mz));
        anRFT(0,0,mz)=tmp;
        tmp=real(anRFT(0,end,mz));
        anRFT(0,end,mz)=tmp;
        if 2*mx==size(anRFT,1)
            tmp=real(anRFT(mx,0,mz));
            anRFT(mx,0,mz)=tmp;
            tmp=real(anRFT(mx,end,mz));
            anRFT(mx,end,mz)=tmp;
        end
    end
else
    error('Unsupported RFT size');
end
outRFT=anRFT; 
