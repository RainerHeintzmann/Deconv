function res=ApplyPhaseModel(vec,afkt)
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
    res=exp(1i*res); % converts the auxilary function back to the all positive object
else  % For data which exists as cell array. E.g. myillu
    savedPhaseInput=res;
    for n=1:size(res,2)
        res{n}=exp(1i*res{n}); % converts the auxilary function back to the all positive object
    end
end

