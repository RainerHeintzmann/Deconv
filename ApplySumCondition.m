% function res=convertVecToIllu(myinput) : Converts a vector of real-space illumination data 
function res=ApplySumCondition(vec, afkt, avoidNeg)
global myillu; % is overwritten!! To save memory
global myim;
global myillu_sumcond;

if nargin < 3
    avoidNeg = 0
end

myillu = afkt(vec); % e.g. convert vector to image using convertVecToObj

asum=0;
sumviews=0;
currentSumCondIdx=1;
for v= 1:numel(myim)  % last pattern will be generated from sum-requirement
    if isempty(myillu_sumcond) || v ~= myillu_sumcond{currentSumCondIdx}
        asum=asum+myillu{v};
        sumviews=sumviews+1;
    else
        tmp=(sumviews+1)-asum;   % The sum is forced to be constrained. 
        if avoidNeg
            tmp(tmp<0)=0;
        end
        myillu{v} = tmp;
        asum=0;sumviews=0;
        currentSumCondIdx=currentSumCondIdx+1;
    end
end
clear asum;

res=myillu;
