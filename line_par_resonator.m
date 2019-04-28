clear variables
close all

FS=12;
model='resonator'

%% Input
Ctuning=1e-9;
nct=100;
Cts=linspace(100e-12,300e-12,nct);    %Range of tuning capacitance
Pnom=5e3;
Rnom=50;
In=sqrt(Pnom/Rnom);
Vn=sqrt(Pnom*Rnom);

%% Computing line parameters
f=50e6;
ro_Cu = 1.72e-8;    %Resistivity of Copper
ro_SS = 6.9e-7;     %Resistivity Stainless Steel 304
ro_Al = 2.62e-8;    %Resistivity Aluminum
ro=1.68e-8;
tand=1e-6;

%230mm out diameter SPINNER
b230=230e-3;        %Outer diameter [m]
a230=100e-3;        %Inner diameter [m]
ro_CuAl = roeq(ro_Cu, ro_Al, a230, b230);
[R230,G230,L230,C230,Z0230] = trasmpar(1,1,tand,ro_CuAl,b230,a230,f);

%150mm out diameter SPINNER
b150=150e-3;
a150=12e-3;
ro_CuAlS = roeq(ro_Cu, ro_Al, a150, b150);
[R150S,G150S,L150S,C150S,Z0150S] = trasmpar(1,1,tand,ro_CuAlS,b150,a150,f);

%150mm out diameter STAINLESS STEEL
b150=150e-3;
a150=12e-3;
ro_CuSS = roeq(ro_Cu, ro_SS, a150, b150);
[R150,G150,L150,C150,Z0150] = trasmpar(1,1,tand,ro_CuSS,b150,a150,f);

%Alumina feedthrough:
b100=100e-3;
a100=12e-3;
epsR1=9.8;
[Rf1,Gf1,Lf1,Cf1,Z0f1] = trasmpar(epsR1,1,tand,ro_CuSS,b100,a100,f);

%Quartz spacers:
epsR2=3.78;
[Rf2,Gf2,Lf2,Cf2,Zf2] = trasmpar(epsR2,1,tand,ro_CuSS,b150,a150,f);


%% Simulink model
set_param('resonator/Arc','commented','on')
yout={};
sps = power_analyze(model,'sort');
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
nfr=5e3;
fr=linspace(35e6,55e6,nfr);
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
[m,n]=max(abs(Velw));
frq=imag(ee)/2/pi;
[MV,IV]=max(m);

decim=10;
cmap=jet(nct/decim);
set(groot,'defaultAxesColorOrder',cmap);

%% Time domain Simulation
Ctuning=Cts(IV);
tsim=1/frq(IV)*1000;
set_param('resonator/Arc','commented','off')
set_param('resonator/Signal Generator','Frequency',num2str(frq(IV)))
simOut=sim(model,'StopTime',num2str(tsim));

%% Plots
figure
ax1=subplot(2,1,1);
grid on
hold on
plot(frq/1e6,m/1e3,'g','linewidth',2)
plot(frq(IV)/1e6,MV/1e3,'o','MarkerFaceColor','w','MarkerSize',10)
ylabel('RMS Voltage [kV]','fontsize',FS);
set(gca,'linewidth',1,'fontsize',FS,'Color','k','XColor','w','YColor','w');
ax2=subplot(2,1,2);
grid on
hold on
plot(frq/1e6,Cts/1e-12,'g','linewidth',2)
plot(frq(IV)/1e6,Cts(IV)/1e-12,'o','MarkerFaceColor','w','MarkerSize',10)
xlabel('Frequency [MHz]','fontsize',FS);
ylabel('Tuning Capacitance [pF]','fontsize',FS);
set(gca,'linewidth',1,'fontsize',FS,'Color','k','XColor','w','YColor','w');
set(gcf,'color','k');
linkaxes([ax1,ax2],'x')

% figure;
% subplot(2,1,1);
% loglog(fr,abs(Zfr(:,1:decim:end)),'linewidth',2);
% grid on;
% colormap(cmap);
% cc=colorbar('ticks',linspace(0,1,nct/decim+1),'ticklabels',[round(Cts(1:decim:end)*1e10)/10, round(Cts(end)*1e10)/10+0.1]);
% cc.Label.String = 'C_{tuning} [nF]';
% xlabel('Frequency [Hz]','fontsize',FS);
% ylabel('Impedance [\Omega]','fontsize',FS);
% title('Total Impedance');
% set(gca,'linewidth',2,'fontsize',FS);
% 
% subplot(2,1,2);
% semilogx(fr,angle(Zfr(:,1:decim:end))*180/pi,'linewidth',2);
% grid on;
% colormap(cmap);
% cc=colorbar('ticks',linspace(0,1,nct/decim+1),'ticklabels',[round(Cts(1:decim:end)*1e10)/10, round(Cts(end)*1e10)/10+0.1]);
% cc.Label.String = 'C_{tuning} [nF]';
% xlabel('Frequency [Hz]','fontsize',FS);
% xlabel('Phase [°]','fontsize',FS);
% set(gca,'linewidth',2,'fontsize',FS);
% 
% % figure;
% % subplot(2,1,1);
% % loglog(fr,abs(Vel),'linewidth',2);
% % grid on;
% % xlabel('Frequency [Hz]','fontsize',FS);
% % ylabel('Impedance [\Omega]','fontsize',FS);
% % set(gca,'linewidth',2,'fontsize',FS);
% % 
% % subplot(2,1,2);
% % semilogx(fr,angle(Vel)*180/pi,'linewidth',2);
% % grid on;
% % xlabel('Frequency [Hz]','fontsize',FS);
% % xlabel('Phase [°]','fontsize',FS);
% % set(gca,'linewidth',2,'fontsize',FS);
% 
% figure;
% semilogx(fr,gain(:,1:decim:end),'linewidth',2);
% grid on;
% colormap(cmap);
% cc=colorbar('ticks',linspace(0,1,nct/decim+1),'ticklabels',[round(Cts(1:decim:end)*1e10)/10, round(Cts(end)*1e10)/10+0.1]);
% cc.Label.String = 'C_{tuning} [nF]';
% xlabel('Frequency [Hz]','fontsize',FS);
% ylabel('Gain [dB]','fontsize',FS);
% title('Ampli to Electrode Voltage Gain')
% set(gca,'linewidth',2,'fontsize',FS);
% 
% set(groot,'defaultAxesColorOrder','remove');
% 
% figure; 
% [C,h1]=contour(CtM,fM,gain,'ShowText','on','LineWidth',4); 
% grid on;
% clabel(C,h1,'fontsize',FS,'labelspacing',500);
% xlabel('C_{tuning}[F]','fontsize',FS);
% ylabel('Frequency [Hz]','fontsize',FS);
% legend('Gain [dB]');
% set(gca,'linewidth',2,'fontsize',FS);
% 
% 
% figure; 
% [C,h1]=contour(CtM,fM,abs(Velw),'ShowText','on','LineWidth',4); 
% clabel(C,h1,'fontsize',FS,'labelspacing',1000);
% xlabel('C_{tuning}[F]','fontsize',FS);
% ylabel('Frequency [Hz]','fontsize',FS);
% legend('Electrode Voltage [V]');
% set(gca,'linewidth',2,'fontsize',FS);
% grid on;
% 
% figure; 
% [C,h1]=contour(CtM,fM,real(Sw),'ShowText','on','LineWidth',4); 
% grid on;
% clabel(C,h1,'fontsize',FS,'labelspacing',1000);
% xlabel('C_{tuning}[F]','fontsize',FS);
% ylabel('Frequency [Hz]','fontsize',FS);
% legend('RF active power [W]');
% set(gca,'linewidth',2,'fontsize',FS);
% 
% figure;
% [C,h1]=contour(CtM,fM,abs(Zfr),0:5:100,'ShowText','on','LineWidth',4); 
% grid on;
% clabel(C,h1,'fontsize',FS,'labelspacing',1000);
% xlabel('C_{tuning}[F]','fontsize',FS);
% ylabel('Frequency [Hz]','fontsize',FS);
% legend('Impdeance [\Omega]');
% set(gca,'linewidth',2,'fontsize',FS);
% 
% figure; 
% [C,h1]=contour(CtM,fM,angle(Zfr)*180/pi,'ShowText','on','LineWidth',4); 
% grid on;
% clabel(C,h1,'fontsize',FS,'labelspacing',1000);
% xlabel('C_{tuning}[F]','fontsize',FS);
% ylabel('Frequency [Hz]','fontsize',FS);
% legend('Phase [°]');
% set(gca,'linewidth',2,'fontsize',FS);







