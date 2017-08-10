function res=ApplyMaskModel(vec,afkt,myRSmask)
global realSpaceMultiplier;
myRSmask=realSpaceMultiplier{myRSmask};
res=afkt(vec);

if ~iscell(res)
    res=myRSmask .* res; % converts the auxilary function back to the all positive object
else  % For data which exists as cell array. E.g. myillu
    for n=1:size(res,2)
        res{n}=myRSmask .* res{n}; % converts the auxilary function back to the all positive object
    end
end

