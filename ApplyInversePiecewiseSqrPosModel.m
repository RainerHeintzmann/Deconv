function res=ApplyInversePiecewiseSqrPosModel(model,afkt)
% model = model +0.0;  % for some CudaMat problem
mask = (model>=1.0);
mask2 = ~mask; 
tmpMask = sqrt(model.*mask + mask2 - 0.75)-0.5;
tmpNotMask = (model - 1.0)./model; % tmp - 1.0
res2 = mask.*tmpMask + mask2.*tmpNotMask;
res = afkt(res2); % the inverse of monotonicPos
