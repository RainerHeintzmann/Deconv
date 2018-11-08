function [myError,residuum]=ResidZero(Recons,myImg,DMask)  % Empty function just returning zero, both for loss and gradient. Useful for testing purposes
myError=0;
residuum=myImg.*0;
