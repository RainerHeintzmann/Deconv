function res=ApplyPhaseModelBwd(thegrad,viewNum,afkt)  % This is obtained using the Wirtinger Derivative of the square model
global savedPhaseInput
% thegrad = 1i*exp(1i*savedPhaseInput{viewNum}).*thegrad;  % To account for the fact that the auxiliary function is what is iterated and the object estimate is the square of it
thegrad = -sin(savedPhaseInput{viewNum}).*real(thegrad) + cos(savedPhaseInput{viewNum}).*imag(thegrad);

res=afkt(thegrad,viewNum);
