function res=ApplyHyperbolicPosModelBwd(thegrad,viewNum,afkt)  % This is obtained using the Wirtinger Derivative of the square model
global savedInputHyperB;
res = thegrad * savedInputHyperB{viewNum};  % This factor was already saved
res=afkt(res,viewNum);
