% res=squeezeIllu(myinput) : converts the cell-arry content in myinput into a squeezed vector not containing the data of the sum condition
function res=squeezeData(myinput)
global my_sumcond;  % 
if isempty(my_sumcond)
    res=cat(4,myinput{:});  % just append everything
else
    asize=size(myinput{1});
    if length(asize) < 3
        asize(length(asize)+1:3)=1;
    end
    res=newim([asize (length(myinput)-length(my_sumcond))]);  % clears the gradient. Defines it first as a dipimage. Later it is converted back to a double vector
    pos=0;currentsumcondidx=1;
    for n=1:size(myinput,2)
        if n ~= my_sumcond{currentsumcondidx}
            res(:,:,:,pos)=myinput{n};
            pos=pos+1;
        else
            currentsumcondidx=currentsumcondidx+1;
        end
    end
end