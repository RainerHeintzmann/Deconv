% function myvec=to3dvec(myvec)  : Forces a size vector to become threedimensional, 3D.
%
function myvec=to3dvec(myvec)
if length(myvec)<3
    myvec(length(myvec)+1:3)=1;
end
