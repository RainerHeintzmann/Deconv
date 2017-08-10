function [myError,residuum,Recons]=ResidRFTLeastSqr(Recons,myImg,DMask)  % a different function is needed here as the sum is computed in RFT space, which is a bit more complicated than in fft space
residuum =(Recons-myImg);

myError = SumAbsSqrRFT(residuum);  % simple least squares, but accounting for complex numbers in residuum

residuum=2*residuum; % to account for the square
