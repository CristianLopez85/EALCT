clear;clc;close all

Hz = 400; tend = 2;
t = 1/Hz : 1/Hz : tend; t = t'; Sf = 100;% time = time(1:400);
% 
amp1 = exp(0.1*t);  amp2 = exp(0.2*t);  amp3 = exp(0.3*t);
% 
Sig1 = amp1.*exp(2*1i*pi*(25*t.^2));                        IF1 = 50*t;
Sig2 = amp2.*exp(2*1i*pi*(25*t));                           IF2 = 25*ones(length(t),1);
Sig3 = amp3.*exp(2*1i*pi*(-22*t.^2 + 95*t + sin(5*t)));     IF3 = -44*t + 95 + 5*cos(5*t);

figure; set(gcf,'Position',[20 100 800 450]);
subplot(221);plot(t,amp1,'b','LineWidth',2);hold on;plot(t,amp2,'r','LineWidth',2);plot(t,amp3,'Color','[0, 0.5, 0]','LineWidth',2);
axis xy; xlim([0,t(end)]);ylim([0.9 1.9]);
legend('Mode 1','Mode 2','Mode 3','Location','northwest')
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');

Sig = Sig1 + Sig2 + Sig3;

% -------------------- STFT ----------------------------
window = 128; Nfrebin = 1024;
[Spec0,f] = STFT(real(Sig(:)),Hz,Nfrebin,window);
subplot(222);imagesc(t,f,abs(Spec0)/max(max(abs(Spec0)))); axis([0 2 0 110]);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'YDir','normal');% set(gca,'FontSize',12);% set(gcf,'Color','w');

%% Estimating the number of components
Sigr = exp(2*pi*1i*(10*t));                       %  IF =  10;
[Spec2,~] = STFT(real(Sigr(:)),Hz,Nfrebin,window);

N = length(t);
n = NaN(1,N);
hop = 8;                                                                   
deltatime = 23;
for i = deltatime : hop: N-(deltatime-1)
    Spec1 = 0*Spec0;
    Spec3 = 0*Spec2;
    pii = i-(deltatime-1); pf = i + (deltatime-1);
    Spec1(:,pii:pf) = Spec0(:,pii:pf);
    H1 = renyi(abs(Spec1));
    Spec3(:,pii:pf) = Spec2(:,pii:pf);
    H2 = renyi(abs(Spec3));
    n(i) = 2^(H1-H2);
end
% 
% figure;plot(floor(n),'bo')
%% ridge extraction and fitting
Num = max(floor(n));    % the number of the components

deltaf = 0.03*Hz; 
beta = 1e-3;
Sig0 = real(Sig)';
if (isreal(Sig0))
    Sig0 = hilbert(Sig0);
end

Nh=151; Nc=71; % threshold=0.03;

[SpecA,Atau,Af] = mALCT(Sig0,Hz,Nh,Nc,Num);
subplot(223);imagesc(Atau,Af,abs(SpecA')/max(max(abs(SpecA')))); axis([0 2 0 110]); axis xy
xlabel('Time (s)'); ylabel('Frequency (Hz)');

cbar = colorbar('horizontal'); % Create horizontal color bar
xlabel(cbar, 'Amplitude [•]'); % Label the color bar
caxis([0 1]);
set(cbar, 'Position', [0.36, 0.46, 0.2, 0.02]);

fidexmult = zeros(Num,length(Sig0)); clearvars Spec Af c bw IFfit extr_Sig findex fidexmult
num = Num;
for i = 1 : Num

    [Spec,Atau,Af] = mALCT(Sig0,Hz,Nh,Nc,num);
    Spec = Spec';
    c = findridges(Spec,deltaf,0.5,num,7);
    bw = 0.01*Hz; % Hz/60 the bandwidth of the TF filter for ridge extraction
    [IFfit,extr_Sig] = Dechirp_filter(Sig0,Hz,bw,Af(c),beta);
    fidexmult(i,:) = c;
    Sig0 = Sig0 - extr_Sig;
    num = num - 1;
end
% figure;plot(Af(fidexmult)');

%% ridge path regrouping

Df = 0.03*Hz;
[findex,interset] = RPRG(fidexmult,Df);

%% Smoothing of ridges

beta1 = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
iniIF = curvesmooth(Af(findex),beta1); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges

subplot(224);
plot(t,IF1,'k','LineWidth',2); hold on; plot(t,IF2,'b','LineWidth',2); plot(t,IF3,'color','[0.75, 0, 0.75]','LineWidth',2);
plot(t,iniIF(3,:),'--c','LineWidth',2);plot(t,iniIF(1,:),'color','0.71 0.71 0.71','LineStyle','--','LineWidth',2);plot(t,iniIF(2,:),'--r','LineWidth',2); 
xlabel('Time (s)'); ylabel('Frequency (Hz)');axis([0 2 0 100]);
legend('T. IF1','T. IF2','T. IF3','E. IF1','E. IF2','E. IF3','Location','north','NumColumns',2)
axis([0 2 0 110]);
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
% figure;,subplot(211);plot(t,iniIF(2,:)-IF1);subplot(212);plot(t,iniIF(1,:)-IF2);
%% run the CT with g for TFC domain
alpha1 = 50;
tt = -8:1/Hz:8;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0
eta = 1;

[tfc1, tfrtic, tcrtic] = sqSTCTm(Sig, eta/length(Sig), 1, h1);

% figure
% imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,201))), 1); axis xy; colormap(1-gray);
% xlabel('frequency (Hz)'); ylabel('chirp rate');
%% TFC ridges
if1 = iniIF(1,:);
if2 = iniIF(2,:);
if3 = iniIF(3,:);
chirpn1 = gradient(if1)/(t(2)-t(1));
chirpn2 = gradient(if2)/(t(2)-t(1));
chirpn3 = gradient(if3)/(t(2)-t(1));

for i = 1:length(Sig)

    [~,bb] = min(abs(tfrtic-if1(i)/Hz));
    line1(i,1) = bb;
    [~,bb1] = min(abs(tcrtic-chirpn1(i)/Hz^2));
    line1(i,2) = bb1-1;

    [~,bb2] = min(abs(tfrtic-if2(i)/Hz));
    line2(i,1) = bb2;
    [~,bb3] = min(abs(tcrtic-chirpn2(i)/Hz^2));
    line2(i,2) = bb3-1;

    [~,bb4] = min(abs(tfrtic-if3(i)/Hz));
    line3(i,1) = bb4;
    [~,bb5] = min(abs(tcrtic-chirpn3(i)/Hz^2));
    line3(i,2) = bb5-1;
end

% reconstruction
scrsz = get(0,'ScreenSize');
rrcon = zeros(3,length(Sig));
scale = alpha1;
stctmat = tfc1;
for char = 1:length(Sig)
    tmp = zeros(3,3);
    tmp(1,1) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn1(char))));
    tmp(1,2) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn2(char))));
    tmp(1,3) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn3(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if3(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn3(char))));

    tmp(2,1) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn1(char))));
    tmp(2,2) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn2(char))));
    tmp(2,3) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn3(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if3(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn3(char))));

    tmp(3,1) = 1/sqrt(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line3(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn1(char))));
    tmp(3,2) = 1/sqrt(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line3(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn2(char))));
    tmp(3,3) = 1/sqrt(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn3(char)))*exp(-pi*(tfrtic(line3(char,1))*Hz-if3(char))^2/(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn3(char))));

    xtmp = [stctmat(line1(char,2),line1(char,1),char); stctmat(line2(char,2),line2(char,1),char); stctmat(line3(char,2),line3(char,1),char)]/Hz;
    rrcon(:,char) = tmp\xtmp;
end
%%
figure; set(gcf,'Position',[20 100 800 450]);
subplot(231)
plot(t,real(Sig1), 'color', 'k'); xlabel('time (s)');hold on
plot(t,real(rrcon(3,:)),'--c');ylim([-2.1 2.1]);
xlabel('Time (s)'); ylabel('Amplitude');%legend('Real','Reconstructed','Orientation','horizontal');

subplot(234)
plot(t, real(Sig1')-real(rrcon(3,:)), 'r'); ylim([-2.1 2.1]); xlabel('Time (s)');ylabel('Amplitude');legend('Error')

subplot(232)
plot(t,real(Sig2), 'color', 'k'); xlabel('Time (s)');hold on;
plot(t,real(rrcon(2,:)),'--c');ylim([-2.1 2.1]);
xlabel('Time (s)'); ylabel('Amplitude');legend('Real','Reconstructed','Orientation','horizontal');

subplot(235)
plot(t, real(Sig2')-real(rrcon(2,:)), 'r'); ylim([-2.1 2.1]); xlabel('Time (s)');ylabel('Amplitude');legend('Error')

subplot(233)
plot(t,real(Sig3), 'color', 'k'); xlabel('Time (s)');hold on;
plot(t,real(rrcon(1,:)),'--c');ylim([-2.1 2.1]);
xlabel('Time (s)'); ylabel('Amplitude');%legend('Real','Reconstructed','Orientation','horizontal');

subplot(236)
plot(t, real(Sig3')-real(rrcon(1,:)), 'r'); ylim([-2.1 2.1]); xlabel('Time (s)');ylabel('Amplitude');legend('Error')

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')