function psf3Db=HalfSidedPSFCycle(psf3D)
    midp=floor(size(psf3D,3)/2);
    if (0) % mod((size(psf3D,3)-midp),2) %Aurelie: not sure how to fix this
        error('This z-size not supported yet.');
    else
        psf3D=psf3D(:,:,midp:end);
        midpn=floor(size(psf3D,3)/2); 
        psf3Db=psf3D;    
        psf3Db(:,:,midpn:end)=psf3D(:,:,0:(size(psf3Db,3)-midpn-1)); %Aurelie: this line may give a bug
        psf3Db(:,:,0:midpn-1)=psf3D(:,:,(size(psf3Db,3)-midpn):(size(psf3Db,3)-midpn)+midpn-1);
    end
    
    