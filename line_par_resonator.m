%impedance reconstruction for IShTAR

clear variables

FS=24;
model='resonator'

Ctuning=1e-9;
nct=100;
Cts=linspace(100e-12,300e-12,nct);
Pnom=5e3;
Rnom=50;
In=sqrt(Pnom/Rnom);
Vn=sqrt(Pnom*Rnom);

%Computing line parameters
f=50e6;
ro=1.68e-8;
tand=1e-6;

%250mm out diameter
b250=230e-3;
a250=100e-3;
d250=0.1;
[R250,G250,L250,C250,Z0250] = trasmpar(1,1,tand,ro,b250,a250,f);
R250=R250*d250;
L250=L250*d250;
C250=C250*d250;

%150mm out diameter
b150=150e-3;
a150=12e-3;
[R150,G150,L150,C150,Z0150] = trasmpar(1,1,tand,ro,b150,a150,f);
% R150=R150*d150;
% L150=L150*d150;
% C150=C150*d150;

%feedthorugh 1:
epsR1=9.8;
df1=0.02;
[Rf1,Gf1,Lf1,Cf1,Z0f1] = trasmpar(epsR1,1,tand,ro,b150,a150,f);
Rf1=Rf1*df1;
Lf1=Lf1*df1;
Cf1=Cf1*df1;

%feedtrough 2:
epsR2=3.78;
df2=0.02;
[Rf2,Gf2,Lf2,Cf2,Zf2] = trasmpar(epsR2,1,tand,ro,b150,a150,f);
Rf2=Rf2*df2;
Lf2=Lf2*df2;
Cf2=Cf2*df2;


yout={};
sps = power_analyze(model,'sort');
%sps_ss = power_analyze('Red_FTEM_4RF_SSM_v1c','ss');
spsd = power_analyze(model,'detailed');
tic
[A,B,C,D,~,~,~,~,~,~,~,~,~] =power_analyze(model);
toc

%recognize right order of output variables:
stramp=find(contains(sps.OutputExpressions,'Vampli'));
strel=find(contains(sps.OutputExpressions,'Vel'));
strct=find(contains(spsd.IndependentStates,'Ctuning'));

[V,De]=eig(A);

ee=diag(De);
fmode=imag(ee)/2/pi;
remode=real(ee);

em_p=ee(fmode>=1e6 & fmode<100e6);
nmp=length(em_p);

V_p=V(:,fmode>=1e6 & fmode<100e6);
II=eye(size(A));
for i=1:nmp
    Ztot(i)=C(stramp,:)/(1i*imag(em_p(i))*II-A)*B;
    
end


% figure;
% surf(abs(V_p),'edgecolor','none');
% grid on;

%frequency response computation
nfr=2e3;
fr=linspace(43e6,48e6,nfr);
w=2*pi*fr;



tic
Zfr=zeros(nfr,nct);
Vel=zeros(nfr,nct);
ee=zeros(nct,1);
At=A;
Bt=B;
it=1;
itmax=100;
dt=10e-6;
G=zeros(nct,1);
for j=1:nct
    
    At(strct,:)=A(strct,:).*(Ctuning./Cts(j));
%     Bt(strct)=B(strct).*(Ctuning./Cts(j));
k=1e-6;
it=1;
st=0;
    while it<=itmax && st==0
    ee(j)=eigs(At+Bt*k*C(stramp,:),1,1i*2*pi*45e6);
    k=abs(k*(1-real(ee(j))*dt));
    k(k>10)=10;
    k(k<1e-6)=1e-6;
    if abs(real(ee(j)))<=1
        st=1;
    end
    it=it+1;
    end
G(j)=k;
    
%     for i=1:nfr
        Zfr(i,j)=C(stramp,:)/(1i*imag(ee(j))*II-At)*Bt;
        Vel(i,j)=C(strel,:)/(1i*imag(ee(j))*II-At)*Bt;
%     end
end
toc

Sw=zeros(nfr,nct);
Vw=zeros(nfr,nct);
Iw=zeros(nfr,nct);
Vw(abs(Zfr)>=50)=Vn;
Sw(abs(Zfr)>=50)=Vn.^2./Zfr(abs(Zfr)>=50);
Sw(abs(Zfr)<50)=In.^2.*Zfr(abs(Zfr)<50);
Vw(abs(Zfr)<50)=Zfr(abs(Zfr)<50).*In;
Iw(abs(Zfr)>=50)=Zfr(abs(Zfr)>=50).\Vn;
Iw(abs(Zfr)<50)=In;
Velw=Vel.*Iw;

[CtM,fM]=meshgrid(Cts,fr);
gain=20*log10(abs(Vel./Zfr));

decim=10;
cmap=jet(nct/decim);
set(groot,'defaultAxesColorOrder',cmap);

figure;
subplot(2,1,1);
loglog(fr,abs(Zfr(:,1:decim:end)),'linewidth',2);
grid on;
colormap(cmap);
cc=colorbar('ticks',linspace(0,1,nct/decim+1),'ticklabels',[round(Cts(1:decim:end)*1e10)/10, round(Cts(end)*1e10)/10+0.1]);
cc.Label.String = 'C_{tuning} [nF]';
xlabel('Frequency [Hz]','fontsize',FS);
ylabel('Impedance [\Omega]','fontsize',FS);
title('Total Impedance');
set(gca,'linewidth',2,'fontsize',FS);

subplot(2,1,2);
semilogx(fr,angle(Zfr(:,1:decim:end))*180/pi,'linewidth',2);
grid on;
colormap(cmap);
cc=colorbar('ticks',linspace(0,1,nct/decim+1),'ticklabels',[round(Cts(1:decim:end)*1e10)/10, round(Cts(end)*1e10)/10+0.1]);
cc.Label.String = 'C_{tuning} [nF]';
xlabel('Frequency [Hz]','fontsize',FS);
xlabel('Phase [°]','fontsize',FS);
set(gca,'linewidth',2,'fontsize',FS);

% figure;
% subplot(2,1,1);
% loglog(fr,abs(Vel),'linewidth',2);
% grid on;
% xlabel('Frequency [Hz]','fontsize',FS);
% ylabel('Impedance [\Omega]','fontsize',FS);
% set(gca,'linewidth',2,'fontsize',FS);
% 
% subplot(2,1,2);
% semilogx(fr,angle(Vel)*180/pi,'linewidth',2);
% grid on;
% xlabel('Frequency [Hz]','fontsize',FS);
% xlabel('Phase [°]','fontsize',FS);
% set(gca,'linewidth',2,'fontsize',FS);

figure;
semilogx(fr,gain(:,1:decim:end),'linewidth',2);
grid on;
colormap(cmap);
cc=colorbar('ticks',linspace(0,1,nct/decim+1),'ticklabels',[round(Cts(1:decim:end)*1e10)/10, round(Cts(end)*1e10)/10+0.1]);
cc.Label.String = 'C_{tuning} [nF]';
xlabel('Frequency [Hz]','fontsize',FS);
ylabel('Gain [dB]','fontsize',FS);
title('Ampli to Electrode Voltage Gain')
set(gca,'linewidth',2,'fontsize',FS);

set(groot,'defaultAxesColorOrder','remove');

figure; 
[C,h1]=contour(CtM,fM,gain,'ShowText','on','LineWidth',4); 
grid on;
clabel(C,h1,'fontsize',FS,'labelspacing',500);
xlabel('C_{tuning}[F]','fontsize',FS);
ylabel('Frequency [Hz]','fontsize',FS);
legend('Gain [dB]');
set(gca,'linewidth',2,'fontsize',FS);


figure; 
[C,h1]=contour(CtM,fM,abs(Velw),'ShowText','on','LineWidth',4); 
clabel(C,h1,'fontsize',FS,'labelspacing',1000);
xlabel('C_{tuning}[F]','fontsize',FS);
ylabel('Frequency [Hz]','fontsize',FS);
legend('Electrode Voltage [V]');
set(gca,'linewidth',2,'fontsize',FS);
grid on;

figure; 
[C,h1]=contour(CtM,fM,real(Sw),'ShowText','on','LineWidth',4); 
grid on;
clabel(C,h1,'fontsize',FS,'labelspacing',1000);
xlabel('C_{tuning}[F]','fontsize',FS);
ylabel('Frequency [Hz]','fontsize',FS);
legend('RF active power [W]');
set(gca,'linewidth',2,'fontsize',FS);

figure;
[C,h1]=contour(CtM,fM,abs(Zfr),0:5:100,'ShowText','on','LineWidth',4); 
grid on;
clabel(C,h1,'fontsize',FS,'labelspacing',1000);
xlabel('C_{tuning}[F]','fontsize',FS);
ylabel('Frequency [Hz]','fontsize',FS);
legend('Impdeance [\Omega]');
set(gca,'linewidth',2,'fontsize',FS);

figure; 
[C,h1]=contour(CtM,fM,angle(Zfr)*180/pi,'ShowText','on','LineWidth',4); 
grid on;
clabel(C,h1,'fontsize',FS,'labelspacing',1000);
xlabel('C_{tuning}[F]','fontsize',FS);
ylabel('Frequency [Hz]','fontsize',FS);
legend('Phase [°]');
set(gca,'linewidth',2,'fontsize',FS);







