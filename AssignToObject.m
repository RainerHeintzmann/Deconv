function thegrad=AssignToObject(res)
global aRecon;
aRecon=res;
if nargout>0
    thegrad=newim(size(res));  % clears the gradient. Defines it first as a dipimage. Later it is converted back to a double vector
end    
