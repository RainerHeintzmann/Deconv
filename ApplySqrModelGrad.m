function res=ApplySqrModelGrad(thegrad,afkt)  % This is obtained using the Wirtinger Derivative of the square model
global savedInput
thegrad = 2*savedInput.*thegrad;  % To account for the fact that the auxilary function is what is iterated and the object estimate is the square of it
% * -10 makes it better but not really good.
res=afkt(thegrad);
