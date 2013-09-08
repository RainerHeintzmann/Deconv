% [intensities,amplitudes,fplanes]=GenSpeckles2(asize,scales,lambda,NA,NumSpeck)  : generates a speckle pattern by random phase values
% asize : image size
% scales : pixel scales (in nm)
% NA : Numerical aperture in air
% NumSpeck : number of speckle patterns to generate
%
% intensities : intensty speckle patterns
% amplitudes : amplitude speckle patterns
% fplanes : Fourier plane speckle patterns

function [intensities,amplitudes,fplanes]=GenSpeckles2(asize,scales,lambda,NA,NumSpeck)  % generates a speckle pattern by random phase values

intensities=cell(NumSpeck,1);
amplitudes=cell(NumSpeck,1);
fplanes=cell(NumSpeck,1);
myrr=rr(asize(1:2),'freq');
myrr2=myrr.*myrr;
for s=1:NumSpeck
    fplane=exp(i*noise(newim(asize(1:2)),'uniform',0,2*pi));
    
    %fplane=ft(fplane);
    %amp=newim(asize(1:2),'scomplex');
    
    sz=1;
    if numel(asize) > 2
        sz=asize(3);
    end
    amp=PropagateFieldTo3D(fplane,scales(3),sz,lambda,scales,NA);
    
    amplitudes{s}=amp;
    intensities{s}=real(amp .* conj(amp));
    fplanes{s}=fplane;
end

