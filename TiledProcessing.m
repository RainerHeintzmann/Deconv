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

Rg=myarglist{5};%fengjiao 25.05.2020
CLEM_GG = 0;
CLEM_IG = 0;
for ri=1:size(Rg,1)
    if iscell(Rg{1})
        Rg = Rg{1};  % in case the parameters are {{ObjParams],[IlluParams],{PSFParams}} use only the ObjParams
    end
    switch (Rg{ri,1})
        case 'CLE_GS'
            CLEM_GG = 1;
            RefX = Rg{ri,2}{2};
            RefY = Rg{ri,2}{3};
        case 'IG'
            CLEM_IG =1 ;
            RefX = Rg{ri,2}{2};
    end
end
if nargin < 5
    aFunction=@GenericDeconvolution
end
if nargin<6
    extraOffset=0;
end
img=data;
clear data
imgsize=size(img{1});
if nargin < 1
    if numel(imgsize)==2 %2D data
        tilesize=[100 100];
    else
        if numel(imgsize)==3 %3D data
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

if any(tilesize>=imgsize)
    fprintf('Warning: Tilesize bigger or equal to image size along a dimension. Removing border in this direction\n')
    aborder(tilesize>=imgsize) = 0;
    fprintf('New border is ');
    fprintf('%d, ',aborder);
    fprintf('\n');
end

if any(aborder>tilesize/2)
    fprintf('Warning: Tile overlap too large along a dimension. Limiting border to 50 percent.\n')
    aborder(aborder>tilesize/2)=floor(tilesize(aborder>tilesize/2)/2);  % limits the border to half the tilesize. floor is important to avoid a negative sized middle part in the cat command.
    fprintf('New border is ');
    fprintf('%d, ',aborder);
    fprintf('\n');
end

tilesize=tilesize-aborder; %for the overlap
tilenum=ceil((imgsize-aborder)./(tilesize)); %number of tiles. Aurelie
aborder(tilenum<=1) = 0;  % Fix! R.H.
tilenum(tilenum <1)=1;
if exist('disableCuda','file')
    disableCuda()
end
resim=newim(imgsize,datatype(img{1})); %initialisation

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
            if iscell(img)
                mytile = cell(1,length(img));
                for n=1:length(img)
                    mytile{n}=img{n}(xs:xe,ys:ye,zs:ze); %Cut a tile. Size is as defined in input, except for the last tile.
                end
            else
                mytile=img(xs:xe,ys:ye,zs:ze); %Cut a tile. Size is as defined in input, except for the last tile.
            end
            if extraOffset~=0
                mytile=mytile+extraOffset;
            end
            if(CLEM_GG)
                myarglist{5}{ri-1,2}{2} = RefX(xs:xe,ys:ye,zs:ze); % fengjiao 11.08.2020
                myarglist{5}{ri-1,2}{3} = RefY(xs:xe,ys:ye,zs:ze); % update the input reference image
            end
            if(CLEM_IG)
                myarglist{5}{ri-1,2}{2} = RefX(xs:xe,ys:ye,zs:ze); %fengjiao 25.05.2020
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

% function res=mycat(dimension,varargin)
% args=varargin(~cellfun('isempty',varargin));
% res = cat(dimension,args{:});


% function res=myramp(sz,dimension,varargin)
% if any(sz<=0)
%     res=[];
% else
%     res=ramp(sz,dimension,varargin{:});
% end
