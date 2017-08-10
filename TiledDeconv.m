% resim=TiledDeconv(tilesize,bordersize,extraOffset,varargin)  : 
% Performs a tiled deconvolution and automatically calculate the required number of tiles and submits them
% individually.
%
% tilesize: 3D size which fits into memory. 
% bordersize: Size of tile-wise overlap.
% extraOffset: Gives an extra offset to each tile.
% varargin: Contains all the other arguments as taken by GenericDeconvolution.

function resim=TiledDeconv(tilesize,bordersize,extraOffset,varargin)
% Goal: perform the deconvolution tile by tile (for the cases where the
% whole dataset is too big for cuda memory), then recombine.
% The tiles should have overlap, then be cut.

resim=TiledProcessing(tilesize,bordersize,varargin{1},{varargin{2:end}},@GenericDeconvolution,extraOffset);

disableCuda();
