function Recons=FwdSubsample(Recons,subSampling)
if any(subSampling~=1)
    Recons=Recons(0:subSampling(1):end,0:subSampling(2):end,0:subSampling(3):end);
end
