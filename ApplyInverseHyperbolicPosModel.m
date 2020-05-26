function res=ApplyInverseHyperbolicPosModel(model,afkt)
b2 = 1.0;
Eps = 1e-9;
model(model<Eps)=Eps; % prevent division by zero
res = afkt((abssqr(model) - b2) / model); % the inverse of monotonicPos
