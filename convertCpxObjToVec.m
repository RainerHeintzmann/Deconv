function out=convertCpxObjToVec(grad)

out=(reshape(double(grad),[prod(size(grad)) 1]));   % this applies to object as well as unpacked (4D) illumination distributions
out=[real(out);imag(out)]; % unpack complex to two reals
