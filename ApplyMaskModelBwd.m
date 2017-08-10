function res=ApplyMaskModelBwd(thegrad,viewNum,afkt,myRSmask)  % This is obtained using the Wirtinger Derivative of the square model
global realSpaceMultiplier;
myRSmask=realSpaceMultiplier{myRSmask};

res = myRSmask.*thegrad;  
res=afkt(res,viewNum);
