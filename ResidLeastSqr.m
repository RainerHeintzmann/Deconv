function [myError,residuum,Recons]=ResidLeastSqr(Recons,myImg,DMask)
residuum =( Recons-myImg);
if ~isempty(DMask)
    residuum(~DMask)=0;
end

myError = real(sum(residuum .* conj(residuum),DMask));  % simple least squares, but accounting for complex numbers in residuum
residuum=2*residuum; % to account for the square
