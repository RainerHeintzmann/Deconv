% res=ConditionalCudaConvert(toConvert,useCuda, forceBin)  : Checks the dataformat wether it is a
% cuda type or not and performs the appropriate converion. Can also handle
% a cell datatype by converting all of the entries.
% toConvert : Array or cell array to convert
% useCuda : 0 : convert to dipimage, 1: convert to cuda
% forceBin: 1 force result to be binary (only for conversion to dipimage)

function  res=ConditionalCudaConvert(toConvert,useCuda,forceBin)
if nargin < 3
    forceBin=0;
end

res=toConvert;

if ~isempty(res) %Aurelie 03.03.2014
    if isa(toConvert,'cell')
        for n=1:numel(toConvert)
            res{n}=ConditionalCudaConvert(res{n},useCuda);
        end
    else
        if isa(toConvert,'cuda') && (useCuda==0)
            res=dip_image_force(res);
        elseif (~isa(toConvert,'cuda')) && (useCuda==1)
            res=cuda(res);
        end
        if forceBin
            res=(res~=0);
        end
    end
end
