% resim=TiledProcessing(tilesize,bordersize,data,myarglist,aFunction,extraOffset)  : Performs a tiled processing of the data
% The first argument is the 3D size which fits into memory. It will then
% automatically calculate the required number of tiles and submit them
% individually to the function.
% aFunction will be called as follows for each tile:
% aFunction(tile,varargin);
% varargin is all other arguments as taken by GenericDeconvolution

function resim=TiledProcessing(tilesize,bordersize,data,myarglist,aFunction,extraOffset)
% Goal: perform the processing tile by tile (for the cases where the
% whole dataset is too big for cuda memory), then recombine.
% The tiles should have overlap, then be cut.
global ringmask;
global DoStop;
ringmask=[];
if nargin < 5
    aFunction=@GenericDeconvolution
end
if nargin<6
    extraOffset=0;
end
img=data;
clear data
imgsize=size(img);
if nargin < 1
    if numel(imgsize)==2 %2D data
        tilesize=[100 100];
    else if numel(imgsize)==3 %3D data
        tilesize=[100 100 100];
        end
    end
end

if length(bordersize) == 1
    if numel(imgsize)==2 %2D data
        aborder=[bordersize bordersize];
    elseif numel(imgsize)==3 %3D data
            aborder=[bordersize bordersize bordersize];
    end
else
    aborder=bordersize;
end

tilesize=tilesize-aborder; %for the overlap
tilenum=ceil((imgsize-aborder)./(tilesize)); %number of tiles. Aurelie
if exist('disableCuda','file')
    disableCuda()
end
resim=newim(img,datatype(img)); %initialisation

tt1=tic;
for tz=1:tilenum(3)
    for ty=1:tilenum(2)
        for tx=1:tilenum(1)    %Perform for each tile
            xs=(tx-1)*tilesize(1); %start index in the image
            xe=min((tx)*tilesize(1)+aborder(1)-1,imgsize(1)-1); %end index in the image. Aurelie
            ys=(ty-1)*tilesize(2);
            ye=min((ty)*tilesize(2)+aborder(2)-1,imgsize(2)-1); %Aurelie
            zs=(tz-1)*tilesize(3);
            ze=min((tz)*tilesize(3)+aborder(3)-1,imgsize(3)-1); %Aurelie
            fprintf('\n\nTiled Processing: Tile [%d,%d,%d] ([%d,%d,%d]), from (%d,%d,%d) to (%d,%d,%d)\n',tx,ty,tz,tilenum,xs,ys,zs,xe,ye,ze);
            mytile=img(xs:xe,ys:ye,zs:ze); %Cut a tile. Size is as defined in input, except for the last tile.
            if extraOffset~=0
                mytile=mytile+extraOffset;
            end
            result=aFunction(mytile,myarglist{:}); %Perform the function on each tile
            disableCuda();
            if isempty(result) || (~isempty(DoStop) && DoStop)
                resim=[];
                return;
            end

            myweights = @(sz,d,left,right) reorient(cat(1,1-cos((ramp(left,1,'freq')+0.5)*pi/2).^2,ramp(sz-left-right,1)*0+1.0,cos((ramp(right,1,'freq')+0.5)*pi/2).^2),d);
            xL = (tx>1) * aborder(1);
            xH = (tx<tilenum(1)) * aborder(1);
            yL = (ty>1) * aborder(2);
            yH = (ty<tilenum(2)) * aborder(2);
            zL = (tz>1) * aborder(3);
            zH = (tz<tilenum(3)) * aborder(3);
            weights = myweights(size(result,1),1,xL,xH).*myweights(size(result,2),2,yL,yH).*myweights(size(result,3),3,zL,zH) ;
            % Now, cut and paste:
            % resim(xs:xe,ys:ye,zs:ze)=result(sxs:sxe,sys:sye,szs:sze); %Aurelie. Also converts the datatype automatically
            resim(xs:xe,ys:ye,zs:ze) = resim(xs:xe,ys:ye,zs:ze) + result * weights; 
            tt=toc(tt1);
            progress=(tx-1+(ty-1)*tilenum(1)+(tz-1)*tilenum(2)*tilenum(1)) / prod(tilenum);
            fprintf('Progress: %g percent, Estimated time to go: %.3g minutes\n',progress*100,(tt/progress)*(1-progress)/60)
        end
    end
end

