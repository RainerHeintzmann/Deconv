function res=ApplyInversePiecewiseSqrPosModel(model,afkt)
model = model +0.0;  % for some CudaMat
mask = (model>=1.0);
mask2 = ~mask;
res2 = model;
if sum(mask) > 0
    tmp=model(mask);
    tmp=sqrt(tmp - 0.75)-0.5;
    res2(mask)=tmp;
end
% tmp = model(mask2)+0.0; % This causes a DipImage Error!!
if sum(mask2) > 0
    tmp = model(mask2);
    tmp1 = (tmp-1.0); % tmp - 1.0
    tmp2 = model(mask2);
    res2(mask2) = tmp1 ./ tmp2; 
end
res = afkt(res2); % the inverse of monotonicPos
