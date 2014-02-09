nump=10000;
t1=noise(nump*(1+0.5*sin(2*pi*xx/8)+0.5*sin(2*pi*yy/4)),'poisson');
t2=noise(nump*(1+0.25*sin(2*pi*xx/8)),'poisson');
t3=noise(nump*(1+0.1*sin(2*pi*xx/8+pi/4)+0.4*sin(2*pi*yy/4+pi/3)),'poisson');

t4=extract(t1,[128 256]);
t5=extract(t3,[128 256]);

cat(3,t4,t5)
figure
c1=frc({t4,t5},[1/8 1/8])

