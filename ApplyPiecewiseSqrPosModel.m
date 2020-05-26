function res=ApplyPiecewiseSqrPosModel(vec,afkt)
global savedInputPiecewiseSqr;
res=afkt(vec);

if ~iscell(res)
    mask = (res>=0);
    mask2 = ~mask;
    savedInputPiecewiseSqr = 1.0 + 2*res;  % this is the parabolic part
    tmp = res(mask2);
    if ~isempty(tmp)
        res2 = 1.0./(1.0-tmp);
        savedInputPiecewiseSqr(mask2) = abssqr(res2);% this is derivative of the hyperbolic part
        res(mask2) = res2;  %  this hyperbola has a value of 1, a slope of 1 and a curvature of 2 at zero X
    end
    tmp2 = res(mask);
    if ~isempty(tmp2)
        tmp = abssqr(tmp2+0.5)+0.75;  % this parabola has a value of 1, a slope of 1  and a curvature of 2 at zero X
        res(mask) = tmp;% The tmp may be needed for cudaMat problems
    end
else  % For data which exists as cell array. E.g. myillu
    savedInputPiecewiseSqr={};
    for n=1:size(res,2)
        resn = res{n} + 0.0;  % CudaMats needs this copy (+0.0) in the DipImage 3 version. Not sure why...
        mask = (resn >= 0);
        mask2 = ~mask;
        tmp = resn(mask2);
        derivFactor = 1.0 + 2*resn;  % this is the parabolic part
        if ~isempty(tmp)
            res2 = 1.0./(1.0-tmp);
            derivFactor(mask2) = abssqr(res2);% this is derivative of the hyperbolic part
            resn(mask2) = res2;  % this parabola has a value of 1, a slope of 1  and a curvature of 2 at zero X
        end
        savedInputPiecewiseSqr{n} = derivFactor;
        tmp2 = resn(mask);
        if ~isempty(tmp2)
            tmp = abssqr(tmp2+0.5)+0.75;  % this parabola has a value of 1, a slope of 1  and a curvature of 2 at zero X
            resn(mask) = tmp;  % The tmp may be needed for cudaMat problems
        end
        res{n} = resn;  %  this hyperbola has a value of 1, a slope of 1 and a curvature of 2 at zero X
    end
end

