% [intensities,amplitudes]=GenSpeckles(asize,relradius,dz,NumSpeck)  : generates a speckle pattern by random phase values
% asize : image size
% relradius : aperture radius
% dz : z distance
% NumSpeck : number of speckle patterns to generate

function [intensities,amplitudes]=GenSpeckles(asize,relradius,dz,NumSpeck)  % generates a speckle pattern by random phase values

intensities=cell(NumSpeck,1);
amplitudes=cell(NumSpeck,1);
for s=1:NumSpeck
    fplane=exp(i*noise(newim(asize(1:2)),'uniform',0,2*pi));
    r2=relradius^2;
    
    %fplane=ft(fplane);
    amp=newim(asize,'scomplex');
    
    sz=1;
    if numel(asize) > 2
        sz=asize(3);
    end
    for n=0:sz-1
        if numel(asize) > 2
            amp(:,:,n)=ift(fplane .* (rr(fplane,'freq')<relradius)*exp(i*2*pi*n*dz*sqrt(r2-rr(fplane,'freq').*rr(fplane,'freq'))));
        else
            amp(:,:)=ift(fplane .* (rr(fplane,'freq')<relradius));
        end
        fprintf('generated speckle pattern %d, plane # %d\n',s,n);
    end
    amplitudes{s}=amp;
    intensities{s}=real(amp .* conj(amp));
end

