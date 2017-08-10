function [myError,residuum]=ResidGaussianWithReadnoise(Recons,myImg,DMask)
global ReadVariance;

XMinusMu=myImg-Recons; % To save some memory
muPlusC=Recons+ReadVariance;

myError =sum(log(muPlusC)+XMinusMu.^2/muPlusC,DMask);  % fast version of Czesar's i-divergence omitting constants

residuum=(Recons .^2 - myImg.^2+Recons+ReadVariance*(1-2*XMinusMu))./((muPlusC).^2);

if ~isempty(DMask)
    residuum(~DMask)=0;  % This assumes always perfect agreement in this area
end
