

obj=readim('Y:\Matlab\images\resolution_coarse.tif')
h=abssqr(ift(rr(obj)<30))
otf=ft(h);otf=otf/max(abs(otf));
pimg=real(ift(ft(obj)*otf));

nimg=noise(pimg/max(pimg)*100,'poisson')

DecNoNoise=GenericDeconvolution(nimg,h,200,'Poisson','LBFGS',{'ForcePos',1},[1 1 1],[0 0 0],[],0);

cat(3,myDeconv,myDeconv1,myDeconv3)

FwdProj=real(ift(ft(DecNoNoise)*otf));


%%
NBoots=20;
Boots=noise(repmat(FwdProj,[1 1 1 NBoots]),'poisson')

for n=1:NBoots
    fprintf('Bootstrapping %d\n',n)
    %GenericDeconvolution(nimg,h,200,'Poisson','RLL',{},[1 1 1],[0 0 0],[],0);
    BootImg=squeeze(Boots(:,:,:,n-1));
    DecBoots{n}=GenericDeconvolution(BootImg,h,200,'Poisson','LBFGS',{'ForcePos',1},[1 1 1],[0 0 0],[],0);
    %DecBoots{n}=GenericDeconvolution(squeeze(Boots(:,:,:,n-1)),h,200,'Poisson','RLL',{},[1 1 1],[0 0 0],[],0);
    % GenericDeconvolution(Boots(:,:,:,n-1),h,200,'LeastSqr',[],{'ForcePos'},[1 1 1],[0 0 0],[],0);
    dipshow(10,DecStarts{n});drawnow();
end
cat(3,DecBoots{:})
%%
NStarts=20;
%StartImg=real(ift(ft(nimg)*otf));
StartImg=0*nimg+mean(nimg);
%DecNoNoise=GenericDeconvolution(nimg,h,200,'Poisson','LBFGS',{'ForcePos',1},[1 1 1],[0 0 0],[],0);
for n=1:NStarts
    fprintf('Start Averaging %d\n',n)
    myStart=noise(StartImg,'gaussian',30);
    myStart=real(ift(ft(myStart)*otf));
    myStart(myStart<1e-5)=1e-5;
    %GenericDeconvolution(nimg,h,200,'Poisson','RLL',{},[1 1 1],[0 0 0],[],0);
    DecStarts{n}=GenericDeconvolution(nimg,h,200,'Poisson','LBFGS',{'StartImg',myStart;'ForcePos',1},[1 1 1],[0 0 0],[],0);
    dipshow(10,DecStarts{n});
    % GenericDeconvolution(Boots(:,:,:,n-1),h,200,'LeastSqr',[],{'ForcePos'},[1 1 1],[0 0 0],[],0);
end
StartMean=mean(cat(3,DecStarts{:}),[],3);
StartVar=var(cat(3,DecStarts{:}),[],3);
cat(3,DecNoNoise,StartMean,StartVar)
%%
BootMean=mean(cat(4,DecBoots{:}),[],4);
BootVar=var(cat(4,DecBoots{:}),[],4);
StartMean=mean(cat(4,DecStarts{:}),[],4);
StartVar=var(cat(4,DecStarts{:}),[],4);

cat(3,DecNoNoise,StartMean,BootMean,BootVar,StartVar)

%%
writeim(nimg,'nimg')
writeim(DecNoNoise,'DeconUnconstrained')  %  5.99615
writeim(FwdProj,'FwdProj')
writeim(obj,'Obj')
writeim(real(ift(otf)),'PSF');
