%FRC   Fourier Ring Correlation
%
% SYNOPSIS:
%  [avgfrc1, avgfrc2] = FRC(in,pixelsizes)
%
% PARAMETERS:
%  in: cell array of input images (at least two required)
%  pixelsizes: vector of pixelsizes (in um)
%
% OUTPUT:
%  avgfsc1: eq(2) from van Heel (sum(F_1F_2^*) / sqrt(sum|F_1|^2  * sum|F_2|^2) )
%  avgfsc2: eq(3) from van Heel (sum(F_1F_2^*) / sum(|F_1||F_2|) )
%
% LITERATURE:
%   M. van Heel, Similarity measures between images, Ultramicroscopy, 21:95-100, 1987. 
%   M. van Heel and M. Schatz, Fourier shell correlation threshold criteria, Journal of Structural Biology, 151:250-262, 2005. 

% (C) Copyright 1999-2010               Pattern Recognition Group
%     All rights reserved               Faculty of Applied Physics
%                                       Delft University of Technology
%                                       Lorentzweg 1
%                                       2628 CJ Delft
%                                       The Netherlands
%
% Bernd Rieger, Feb 2010
% May 2011, bug fix in normalization
% Rainer Heintzmann, Nov 2013

function [avgfrc1, avgfrc2, binpos] = frc(in,pixelsizes)

% Avoid being in menu
if nargin == 1 & ischar(in) & strcmp(in,'DIP_GetParamList')
   avgfrc1 = struct('menu','none');
   return
end

%if ndims(in)~=3; error('dpr requires at least two images combined in a 3D image');end
if ~iscell(in) ||numel(in) < 2 ; error('frc requires at least two images combined in a cell array');end
% if size(in{1},in{2}) ~=size(in,2); error('Image must be square to perform angular average easily.');end

innerRadius=0;
%sz = size(in);
%N = sz(3);
N=numel(in);

resfrc1 = 0; % newim(max(size(in{1}))/2+1,N-1);
resfrc2 = 0; % newim(max(size(in{1}))/2+1,N-1);
fpixelsizes=1./pixelsizes ./ size(in{1});
dq=min(fpixelsizes)*2;
for ii=0:N-2
   f1 = ft(in{ii+1});
   f2 = ft(in{ii+2});
   m1 = abs(f1);
   m2 = abs(f2);
   [tmp,binpos] = radialmean2( f1*conj(f2),fpixelsizes,dq,innerRadius);
   tmp=real(tmp);
   resfrc1 = resfrc1 + tmp / sqrt( radialmean2(m1^2,fpixelsizes,dq,innerRadius) .* radialmean2(m2^2,fpixelsizes,dq,innerRadius) );
   resfrc2 = resfrc2 + tmp / radialmean2(m1*m2,fpixelsizes,dq,innerRadius);
end
avgfrc1 = resfrc1/(N-1);
avgfrc2 = resfrc2/(N-1);

plot(binpos,avgfrc1);
xlabel('Spatial Frequency [µm^-1]');
ylabel('Correlation Coefficient');
