function res=ApplyMaskModelGrad(thegrad,afkt,myRSmask)  % This is obtained using the Wirtinger Derivative of the square model
global realSpaceMultiplier;
myRSmask=realSpaceMultiplier{myRSmask};

thegrad = myRSmask.*thegrad;  
res=afkt(thegrad);
