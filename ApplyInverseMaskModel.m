function res=ApplyInverseMaskModel(model,afkt,myRSmask)
global realSpaceMultiplier;
myRSmask=realSpaceMultiplier{myRSmask};

res = model ./ myRSmask;  
res=afkt(res);
