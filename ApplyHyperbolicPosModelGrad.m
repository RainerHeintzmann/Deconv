function res=ApplyHyperbolicPosModelGrad(thegrad,afkt)  % This is obtained using the Wirtinger Derivative of the square model
global savedInputHyperB
thegrad = 2*savedInputHyperB.*thegrad;  % To account for the fact that the auxilary function is what is iterated and the object estimate is the square of it
% * -10 makes it better but not really good.
res=afkt(thegrad);
