function res=ApplyInversePiecewiseSqrPosModel(model,afkt)
mask=model>=1.0;
mask2 = ~mask;
res2=model * 0.0;
res2(mask)=sqrt(model(mask) - 0.75)-0.5;
res2(mask2)= (model(mask2)-1.0)./model(mask2);
res = afkt(res2); % the inverse of monotonicPos
