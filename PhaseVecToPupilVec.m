function res=PhaseVecToPupilVec(vec,afkt)
global AmpFactor; 

% global mytmp; %Variable modified in AmpMaskToIllu. Aurelie 03.03.2014
% mytmp={}; %clear in case of further occurances. Aurelie 03.03.2014
if nargin > 1
   res=afkt(vec);
else
   res=vec;
end
global savedPhaseInput;
if ~iscell(res)
    savedPhaseInput={res};
    resR=AmpFactor*cos(res); % converts the auxilary function back to the all positive object
    resI=AmpFactor*sin(res); % converts the auxilary function back to the all positive object
    res=[resR;resI]; % pack into the format of a complex pupil vector
else  % For data which exists as cell array. E.g. myillu
    savedPhaseInput=res;
    for n=1:size(res,2)
        resR=AmpFactor*cos(res{n}); % converts the auxilary function back to the all positive object
        resI=AmpFactor*sin(res{n}); % converts the auxilary function back to the all positive object
        res{n}=[resR;resI]; % pack into the format of a complex pupil vector
    end
end
