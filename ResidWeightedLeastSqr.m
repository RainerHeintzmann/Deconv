function [myError,residuum,Recons]=ResidWeightedLeastSqr(Recons,myImg,DMask)
global DeconvVariances;

residuum =( Recons-myImg);
if ~isempty(DMask)
    residuum(~DMask)=0;
end

if ~isempty(DeconvVariances)
    myError = real(sum(residuum .* conj(residuum) ./ DeconvVariances,DMask));  % simple least squares, but accounting for complex numbers in residuum
    residuum=2*residuum ./ DeconvVariances;
else
    error('If using the WeightedLeastSqr method, you need ot provide a non-empty variance array. (second to last paramter of the call);');
end
