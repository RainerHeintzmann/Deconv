function [myError,residuum]=ResidPoisson(Recons,myImg,DMask)
% The first one below is numerically very unstable
% myError =sum(Recons-myImg .* log(Recons));  % fast version of Czesar's i-divergence omitting constants
Recons=Recons+0;  % TO PREVENT MODIFICATION OF Recons as an argument.
eps=1e-8;  % 1e-7 . To avoid to take log(0) in the Sterling approximation. If eps is too small, the ratio can generate quite high values in the ratio
amask=(Recons<eps);
Recons(amask)=eps; % Modifies the data!
clear amask;
ratio = myImg ./ Recons; 
epsR=1e-7;
amask=(ratio<epsR);
ratio(amask)=epsR;
clear amask;
%ratio(:,:,23)
if ~isempty(DMask)
    ratio(~DMask)=1;  % This assumes always perfect agreement in this area
end

Recons=Recons-myImg; % To save some memory
%myError =sum((Recons-myImg)+myImg .* log(ratio));  % fast version of Czesar's i-divergence omitting constants
myError =sum(Recons+myImg .* log(ratio),DMask);  % The "+" is due to taking the ratio rather than Recons directly
%myError = 0;
%myError =sum(Recons+myImg .* log(ratio));  % 
%myError =sum(Recons-myImg .* log(Recons));  % fast version of Czesar's i-divergence omitting constants

residuum=1-ratio;
