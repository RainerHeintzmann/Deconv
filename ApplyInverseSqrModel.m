function res=ApplyInverseSqrModel(model,afkt)
model = model +0.0;  % needed for CudaMat to make a copy
mask = model<0;
model(mask)=0; % prevent negative values
res=afkt(sqrt(model));
