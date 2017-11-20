function startVec=PreProcessStartVec(startVec,myim,otfrep)
global aRecon;
    if (1)
        svecsize=NewIlluSize; % as defined above.  Not size(myim{1}) is the object lives in a different space
        if RegObj(7,1)  % this means a complex object shall be reconstructed. Thus the start vector also needs to be complex
            startVec=newim(svecsize,'scomplex');
        else
            startVec=newim(svecsize);
        end
        if RegObj(9,1) && (mymean < 1e-2)   % Force Positive. In this case one needs an offset, even though the model (e.g. in the amplite world) does not support it.
            startVec=startVec+1e-3;
            % startVec=NumViews*(double(repmat(1e-3,[1 prod(svecsize)])))';  % For now, just start with a guess of one, so the amplitude algorithm does not screw up
        else
            startVec=startVec+mymean-RegObj(12,1);  % account for the background value in the object estimate
            % startVec=NumViews*(double(repmat(mymean,[1 prod(svecsize)])))';
        end
        aRecon=startVec;
    else
        global para;
        startVec=NumViews*para.res_object_resampled * mymean / mean(para.res_object_resampled)';   % To test it on SIM simulations
    end    
    
    if (~RegObj(7,1))  % Even though the data may be complex, the model may be forced to be real.
        startVec=abs(startVec); % sqrt will be applied later if needed. This is only to force it to be a real valued quantity even if the data is complex
    end
    if (RegObj(8,1)) % IntensityData
        for v=1:length(myim)
            if ~isreal(myim{1})
                error('For intensity data as specified in the flag, all images need to be real valued!');
            end
        end
        % startVec=startVec*exp(i*pi/4);  % Maybe random phases should be supplied as a start?
        if RegObj(7,1)  % 'complex'
            startVec=sqrt(startVec)+1e-8*1i;  % Maybe random phases should be supplied as a start?
        else
            startVec=sqrt(startVec);  % Maybe random phases should be supplied as a start?
        end
        if ~RegObj(21,1)  % not FTData
            startVec=startVec/sqrt(sum(otfrep{1}));
        end
    end

    