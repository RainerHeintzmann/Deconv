function out=convertCpxObjToRVec(grad)

out=(reshape(double(real(grad)),[prod(size(grad)) 1]));   % this applies to object as well as unpacked (4D) illumination distributions
