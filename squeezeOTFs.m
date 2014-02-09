% illuvec=squeezeIllu() : converts the content in the global myillu variable into a squeezed vector not containing the data of the sum condition
function illuvec=squeezeIllu()
global myillu;
global myillu_sumcond;
asize=size(myillu{1});
if length(asize) < 3
    asize(length(asize)+1:3)=1;
end
illuvec=newim([asize (length(myillu)-length(myillu_sumcond))]);  % clears the gradient. Defines it first as a dipimage. Later it is converted back to a double vector
pos=0;currentsumcondidx=1;
for n=1:size(myillu,2)
    if n ~= myillu_sumcond{currentsumcondidx}
        illuvec(:,:,:,pos)=myillu{n};
        pos=pos+1;
    else
        currentsumcondidx=currentsumcondidx+1;
    end
end
