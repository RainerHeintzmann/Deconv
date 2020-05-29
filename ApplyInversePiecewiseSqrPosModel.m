function res=ApplyInversePiecewiseSqrPosModel(model,afkt)
model = model +0.0;  % for some CudaMat problem
mask = (model>=1.0)+0.0;
mask2 = ~mask; % for some CudaMat problem
sum1 = sum(mask);
sum2 = sum(~mask);
res2 = model;
if sum1 > 0
    tmp=model(mask);
    tmp2=sqrt(tmp - 0.75)-0.5;
    clear tmp;
    res2(mask)=tmp2;
end
% tmp = model(mask2)+0.0; % This causes a DipImage Error!!
try
    res2 = res2+0.0;
catch
    % do nothing. Just to catch the error
end

if sum2 > 0
    tmp = model(mask2);
    tmp1 = (tmp-1.0); % tmp - 1.0
    tmp2 = model(mask2);
    res2(mask2) = tmp1 ./ tmp2; 
end
res = afkt(res2); % the inverse of monotonicPos
