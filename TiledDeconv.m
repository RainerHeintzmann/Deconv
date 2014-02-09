% resim=TiledDeconv(tilesize,varargin)  : Performs a tiled deconvolution.
% The firs argument is the 3D size which fits into memory. It will then
% automatically calculate the required number of tiles and submit them
% individually.
% varargin is all other arguments as taken by GenericDeconvolution

function resim=TiledDeconv(tilesize,bordersize,extraOffset,varargin)
% Goal: perform the deconvolution tile by tile (for the cases where the
% whole dataset is too big for cuda memory), then recombine.
% The tiles should have overlap, then be cut.

resim=TiledProcessing(tilesize,bordersize,varargin{1},{varargin{2:end}},@GenericDeconvolution,extraOffset);

disableCuda();
