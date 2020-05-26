function res=ApplyPiecewiseSqrPosModelBwd(thegrad,viewNum,afkt)  % This is obtained using the Wirtinger Derivative of the square model
global savedInputPiecewiseSqr;
res = thegrad * savedInputPiecewiseSqr{viewNum};  % This factor was already saved
res=afkt(res,viewNum);
