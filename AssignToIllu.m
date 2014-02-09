function thegrad=AssignToIllu(res)
global myillu;
global myillu_sumcond;
myillu=res;

asize=size(myillu{1});
if length(asize) < 3
    asize(length(asize)+1:3)=1;
end

if nargout > 0
    thegrad=newim([asize (length(myillu)-length(myillu_sumcond))]);  % clears the gradient. Defines it first as a dipimage. Later it is converted back to a double vector
end
