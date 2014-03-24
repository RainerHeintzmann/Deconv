%% Script that performs all the test
addpath('Y:\MATLAB\Toolboxes\Deconv\GradientTests');

global todisplay

%Illu gradient with mask
    %Without positivity constraint
    
    %No cuda
    ErrIlluMaskNoCudaNoPos=PerformTestGradientGen('myType','illu','useMask',1,'NumIm',4)
    %  10.03.2013: 2.6705e-04
    Display_illuMaskNoCudaNoPos=todisplay;
    %Cuda
    ErrIlluMaskCudaNoPos=PerformTestGradientGen('myType','illu','useMask',1,'NumIm',4,'useCuda',1)
    % 10.03.2013: 1.2526
    Display_illuMaskCudaNoPos=todisplay;

    %With positivity constraint
    
    %cuda
    ErrIlluMaskCudaPos=PerformTestGradientGen('myType','illu','useMask',1,'NumIm',4,'useCuda',1,'MyReg',{'ForcePos',[]})
    % 10.03.2013: 0.8318
    Display_illuMaskCudaPos=todisplay;
    %no cuda
    ErrIlluMaskNoCudaPos=PerformTestGradientGen('myType','illu','useMask',1,'NumIm',4,'useCuda',0,'MyReg',{'ForcePos',[]})
    % 10.03.2013: 2.5196e-04
    Display_illuMaskNoCudaPos=todisplay;

%% The rest has not been tested yet!
%Illu gradient without mask
    %Without positivity constraint
    
    %No cuda
    ErrIlluNoMaskNoCudaNoPos=PerformTestGradientGen('myType','illu','NumIm',2,'sX',10,'sY',10,'useCuda',0);
    Display_illuMaskNoCudaPos=todisplay;
    %Cuda
    
    %With positivity constraint
    
    %No cuda
    ErrIlluNoMaskNoCudaPos=PerformTestGradientGen('myType','illu','NumIm',2,'sX',10,'sY',10,'useCuda',0,'MyReg',{'ForcePos',[]});
    Display_illuMaskNoCudaPos=todisplay;
    %Cuda
