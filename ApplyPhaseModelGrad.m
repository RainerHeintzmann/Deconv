function res=ApplyPhaseModelGrad(thegrad,afkt)  % This is obtained using the Wirtinger Derivative of the square model
global savedPhaseInput
thegrad = exp(i*savedPhaseInput).*thegrad;  % To account for the fact that the auxilary function is what is iterated and the object estimate is the square of it

res=afkt(thegrad);
