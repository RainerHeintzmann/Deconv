function [normImg, myNorm] = FwdBlindIntensity(img,measInt)
global ComplexData
if ComplexData == 1
    myNorm = sum(abs(img).^2);
else
    myNorm = sum(img); 
end
normImg = img*measInt/myNorm;