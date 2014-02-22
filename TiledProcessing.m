% resim=TiledProcessing(tilesize,data,varargin,aFunction)  : Performs a tiled processing of the data
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
else if numel(imgsize)==3 %3D data
    aborder=[bordersize bordersize bordersize];
    end
end
else
    aborder=bordersize;
end

tilesize=tilesize-aborder; %for the overlap
tilenum=ceil((imgsize-aborder)./(tilesize)); %number of tiles. Aurelie
resim=newim(img,datatype(img)); %initialisation

for tz=1:tilenum(3)
    for ty=1:tilenum(2)
        for tx=1:tilenum(1)    %Perform for each tile
            xs=(tx-1)*tilesize(1); %start index in the image
            xe=min((tx)*tilesize(1)+aborder(1)-1,imgsize(1)-1); %end index in the image. Aurelie
            ys=(ty-1)*tilesize(2);
            ye=min((ty)*tilesize(2)+aborder(2)-1,imgsize(2)-1); %Aurelie
            zs=(tz-1)*tilesize(3);
            ze=min((tz)*tilesize(3)+aborder(3)-1,imgsize(3)-1); %Aurelie
            fprintf('Tiled Processing: Tile [%d,%d,%d] ([%d,%d,%d]), from (%d,%d,%d) to (%d,%d,%d)\n',tx,ty,tz,tilenum,xs,ys,zs,xe,ye,ze);
            mytile=img(xs:xe,ys:ye,zs:ze); %Cut a tile. Size is as defined in input, except for the last tile.
            if extraOffset~=0
                mytile=mytile+extraOffset;
            end
            result=aFunction(mytile,myarglist{:}); %Perform the function on each tile
            sxs=0; %start index of the tile
            sys=0;
            szs=0;
            sxe=size(result,1)-1; %end index of the tile. Aurelie
            sye=size(result,2)-1; %Aurelie
            sze=size(result,3)-1; %Aurelie
            if tx>1
                xs=xs+floor(aborder(1))/2; %start index to write on final image
                sxs=floor(aborder(1))/2; %cut away border (start index to read result)
            end
            if ty>1
                ys=ys+floor(aborder(2))/2;
                sys=floor(aborder(2))/2;
            end
            if tz>1
                zs=zs+floor(aborder(3))/2;
                szs=floor(aborder(3))/2;
            end
            if tx<tilenum(1)
                xe=xe-floor(aborder(1))/2; %end index to write on final image
                sxe=sxe-floor(aborder(1))/2; %cut away border (end index to read result). Aurelie
            end
            if ty<tilenum(2)
                ye=ye-floor(aborder(2))/2;
                sye=sye-floor(aborder(2))/2; %Aurelie
            end
            if tz<tilenum(3)
                ze=ze-floor(aborder(3))/2;
                sze=sze-floor(aborder(3))/2; %Aurelie
            end
            % Now, cut and paste:
            resim(xs:xe,ys:ye,zs:ze)=result(sxs:sxe,sys:sye,szs:sze); %Aurelie. Also converts the datatype automatically
        end
    end
end
