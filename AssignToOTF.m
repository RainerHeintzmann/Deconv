function thegrad=AssignToOTF(res) % input will be an otf
global otfrep;
global aRecon;
otfrep=res;
DataSize=size(aRecon{1});  % needs to be for the PSF domain
if nargout > 0
    thegrad=newim([DataSize length(otfrep)]);  % clears the gradient. Defines it first as a dipimage. Later it is converted back to a double vector
end
