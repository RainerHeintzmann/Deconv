function CalcBwdOTF(paramVec) % calculate the Wiener-Butterworth OTF according to the Paper from the shroff group: https://www.biorxiv.org/content/10.1101/647370v1.supplementary-material
global BwdOTF;
global otfrep;
% find the CutOff directly from the OTF.
if length(paramVec) < 1 || paramVec(1) == 0
    order = 7;
else
    order = paramVec(1);
end
if length(paramVec) < 2 || paramVec(2) == 0
    alpha = 0.001;
else
    alpha = paramVec(2);
end
if length(paramVec) < 3 || paramVec(3) == 0
    beta = 0.001;
else
    beta = paramVec(3);
end
for n=1:length(otfrep)
    OTF = otfrep{n};
    aOTF = abs(OTF);
    clear('OTF');
    ZeroTransfer = max(aOTF);
    aOTF = aOTF / ZeroTransfer;
    RelBorder = 1e-3;
    if (0)
        ProjX = sum(aOTF,[],[2,3]); ProjX = ProjX / max(ProjX);ProjY = sum(aOTF,[],[1,3]); ProjY = ProjY / max(ProjY);ProjZ = sum(aOTF,[],[1,2]); ProjZ = ProjZ / max(ProjZ);
        XMask = ProjX < RelBorder;YMask = ProjY < RelBorder;ZMask = ProjZ < RelBorder;
        XPos = find(XMask,1);YPos = find(YMask,1);ZPos = find(ZMask,1);
    else
        OTFMask = aOTF > RelBorder;
        ProjX = sum(OTFMask,[],[2,3]); ProjX = ProjX / max(abs(ProjX));ProjY = sum(OTFMask,[],[1,3]); ProjY = ProjY / max(abs(ProjY));ProjZ = sum(OTFMask,[],[1,2]); ProjZ = ProjZ / max(abs(ProjZ));
        if isa(ProjX,'cuda')
            ProjX=double_force(ProjX);
            ProjY=double_force(ProjY);
            ProjZ=double_force(ProjZ);
        end
        XPos = find(~ProjX,1);YPos = find(~ProjY,1);ZPos = find(~ProjZ,1);
    end
    CutOff = [XPos,YPos,ZPos];  % Abbe Radius.  AbbeMaxRadiusFromPSF(ImageParam,PSFParam);
    CutOffRed = CutOff*0.5;
    OTFVal = (aOTF(round(CutOffRed(1)),0,0) + aOTF(0,round(CutOffRed(2)),0) + aOTF(0,0,round(CutOffRed(3)))) / 3.0;  % This is the CutOffGain
    fprintf('Wiener-Butterworth: orders: %d, alpha: %g, beta: %g, cutOffGain: %g\n',order,alpha,beta,double(OTFVal));
    BWiener = OTFVal;
    % BWiener = OTFVal/(abssqr(OTFVal)+alpha);
    % eps2 = abssqr(BWiener/beta)-1.0;  % eps^2, see equation 27 in https://www.biorxiv.org/content/10.1101/647370v1.supplementary-material
    eps2 = BWiener/abssqr(beta)-1.0;  % eps^2, see equation 27 in https://www.biorxiv.org/content/10.1101/647370v1.supplementary-material
    Wiener = (aOTF)/(abssqr(aOTF)+alpha)*(1+alpha);  % the "conj" has been omitted since it is included in the Bwd projection
    sOTF = (size(aOTF)-[0,1,0]) .* [1,2,1];
    Butterworth = fft2rft(ifftshift(1.0/sqrt(1.0+eps2*rrscale(sOTF,1.0./CutOff)^(2.0*order))));
    BwdOTF{n} = Wiener .* Butterworth;
    % figure; plot(abs(BwdOTF{1}(:,0,0))); hold on; plot(abs(Wiener(:,0,0)),'g');plot(abs(Butterworth(:,0,0)),'r');
    BwdOTF{n} = ZeroTransfer * BwdOTF{n};
end
