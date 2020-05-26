function res=ApplyPiecewiseSqrPosModel(vec,afkt)
global savedInputPiecewiseSqr;
res=afkt(vec);

if ~iscell(res)
    mask = res>=0;
    mask2 = ~mask;
    savedInputPiecewiseSqr = 1.0 + 2*res;  % this is the parabolic part
    res2 = 1.0./(1.0-res(mask2));
    savedInputPiecewiseSqr(mask2) = abssqr(res2);% this is derivative of the hyperbolic part
    res(mask2) = res2;  %  this hyperbola has a value of 1, a slope of 1 and a curvature of 2 at zero X
    res(mask) = abssqr(res(mask)+0.5)+0.75;  % this parabola has a value of 1, a slope of 1  and a curvature of 2 at zero X
else  % For data which exists as cell array. E.g. myillu
    savedInputPiecewiseSqr={};
    for n=1:size(res,2)
        resn = res{n};
        mask = res{n}>=0;
        mask2 = ~mask;
        res2 = 1.0./(1.0-resn(mask2));
        derivFactor = 1.0 + 2*resn;  % this is the parabolic part
        derivFactor(mask2) = abssqr(res2);% this is derivative of the hyperbolic part
        savedInputPiecewiseSqr{n} = derivFactor;
        resn(mask) = abssqr(resn(mask)+0.5)+0.75;  % this parabola has a value of 1, a slope of 1  and a curvature of 2 at zero X
        resn(mask2) = res2;  % this parabola has a value of 1, a slope of 1  and a curvature of 2 at zero X
        res{n} = resn;  %  this hyperbola has a value of 1, a slope of 1 and a curvature of 2 at zero X
    end
end

