function res=ApplyInverseSqrModel(model,afkt)
model(model<0)=0; % prevent negative values
res=afkt(sqrt(model));
