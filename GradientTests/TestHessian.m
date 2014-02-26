aRecon=gaussf((abs(xx+yy)< 10) * (abs(yy) < 50));
aRecon(20,20)=1;

s1=1/1.6;s2=0.1727;
H11=dip_convolve1d(aRecon,s1*[1 -2 1],0,1); % second X derivative
H22=dip_convolve1d(aRecon,s1*[1 -2 1],1,1); % second Y derivative
H2a=dip_convolve1d(aRecon,[-1 0 1],0,1); %
H12=dip_convolve1d(H2a,s2*[-1 0 1],1,1); % mixed derivative XY
fprintf('H11+H22: %g, 2*H12 = %g, Hges: %g\n',sum(H11 .*H11 + H22.*H22), sum(2*H12.*H12),sum(H11 .*H11 + H22.*H22+2*H12.*H12))  % Does it need the weights of the filters?
H=hessian(aRecon);
fprintf('H11+H22: %g, 2*H12 = %g, Hges: %g\n',sum(H{1,1}*H{1,1}+H{2,2}*H{2,2}), sum(2*H{1,2}*H{1,2}),sum(H{1,1}*H{1,1}+H{2,2}*H{2,2}+2*H{1,2}*H{1,2}))  % Does it need the weights of the filters?

aRecon=rotation(aRecon*1.0,pi/8);
H11=dip_convolve1d(aRecon,s1*[1 -2 1],0,1); % second X derivative
H22=dip_convolve1d(aRecon,s1*[1 -2 1],1,1); % second Y derivative
H2a=dip_convolve1d(aRecon,[-1 0 1],0,1); %
H12=dip_convolve1d(H2a,s2*[-1 0 1],1,1); % mixed derivative XY
fprintf('H11+H22: %g, 2*H12 = %g, Hges: %g\n',sum(H11 .*H11 + H22.*H22), sum(2*H12.*H12),sum(H11 .*H11 + H22.*H22+2*H12.*H12))  % Does it need the weights of the filters?

H=hessian(aRecon);
fprintf('H11+H22: %g, 2*H12 = %g, Hges: %g\n',sum(H{1,1}*H{1,1}+H{2,2}*H{2,2}), sum(2*H{1,2}*H{1,2}),sum(H{1,1}*H{1,1}+H{2,2}*H{2,2}+2*H{1,2}*H{1,2}))  % Does it need the weights of the filters?


