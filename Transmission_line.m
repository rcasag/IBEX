
c=3e8;
f=50e6;
w=2*pi*f;
B=w/c;
lambda=c/f;
t=linspace(0,100,1000);
y=1/(2i)*(exp(1i*w*t)-exp(-1i*w*t));
plot(t,y)