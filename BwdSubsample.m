function SubSampled=BwdSubsample(Recons,subSampling)
if any(subSampling~=1)
    SubSampled=newim(size(Recons).*subSampling);
    SubSampled(0:subSampling(1):end,0:subSampling(2):end,0:subSampling(3):end)=Recons;
else
    SubSampled=Recons;
end
