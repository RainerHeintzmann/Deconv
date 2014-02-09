function res=ApplySqrModel(vec,afkt)
res=afkt(vec);

global savedInput;
if ~iscell(res)
    savedInput=res;
    res=abssqr(res); % converts the auxilary function back to the all positive object
else  % For data which exists as cell array. E.g. myillu
    savedInput=squeezeData(res);
    for n=1:size(res,2)
        res{n}=abssqr(res{n}); % converts the auxilary function back to the all positive object
    end
end

