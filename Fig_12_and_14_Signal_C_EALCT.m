clc; clear; close all
Hz = 100;
t = 0:1/Hz:3; t = t';

am1 = ones(1,length(t)); am2 = 0.9*am1';

figure; set(gcf,'Position',[20 100 800 450]);
subplot(221);plot(t,am2,'color','0.25 0.80 0.54','LineWidth',2);hold on;plot(t,am1,'color','0.5 0 1','LineWidth',2);
axis xy; xlim([0,t(end)]);ylim([0.6 1.4]);
legend('Mode 1','Mode 2','Location','northwest')
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');

Sig1 = exp(2*pi*1i*(30*t + 3*sin(3*t)));            IF2 = 30 + 9*cos(3*t);
Sig2 = 0.9*exp(2*pi*1i*(30*t - 3*sin(3*t) + 3/(2*pi)));            IF1 = 30 - 9*cos(3*t);

initstate(2) ;
Sig = Sig1 + Sig2;
Sig = awgn(Sig,10);

% figure; set(gcf,'Position',[20 100 800 450]); subplot(221);plot(t,real(Sig),'k');
% xlabel('Time (s)'); ylabel('Amplitude');

% figure
% set(gcf,'Position',[20 100 320 250]);% set(gcf,'Color','w');
% plot(t,Sig,'linewidth',1);
% xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
% ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');
% set(gca,'FontSize',12);% set(gca,'linewidth',1);

%%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
window = 128/2; Nfrebin = 1024;
[Spec0,f] = STFT(real(Sig(:)),Hz,Nfrebin,window);
subplot(222);imagesc(t,f,abs(Spec0)/max(max(abs(Spec0)))); axis([0 3 10 50]);

xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'YDir','normal');% set(gca,'FontSize',12);% set(gcf,'Color','w');
%% Estimating the number of components
E = sum(abs(Spec0)); mu = min(E)/max(E);       % ---- Added by CL, from Miao
D_E = 0.5*mu*max(max(abs(Spec0)));                 % ---- Added by CL, from Miao

Sigr = exp(2*pi*1i*(10*t));                % IF2 =  10;
[Spec2,~] = STFT(real(Sigr(:)),Hz,Nfrebin,window);

N = length(t);
n = NaN(1,N);
deltatime = 23;
for i = deltatime : N-(deltatime-1)
    Spec1 = 0*Spec0;
    Spec3 = 0*Spec2;
    pii = i-(deltatime-1); pf = i + (deltatime-1);
    Spec1(:,pii:pf) = Spec0(:,pii:pf);
    Spec1(abs(Spec1)<D_E) = 0;
    H1 = renyi(abs(Spec1));
    Spec3(:,pii:pf) = Spec2(:,pii:pf);
    H2 = renyi(abs(Spec3));
    Spec3(abs(Spec3)<D_E) = 0;
    n(i) = 2^(H1-H2);
end

figure;plot(floor(n),'bo')
%% ridge extraction and fitting
Num = max(floor(n));    % the number of the components

deltaf = 0.05*Hz;
beta = 1e-3;
Sig0 = real(Sig)';
if (isreal(Sig0))
    Sig0 = hilbert(Sig0);
end

Nh=71; Nc=31; % threshold=0.03;

[SpecA,Atau,Af] = mALCT(Sig0,Hz,Nh,Nc,Num);
subplot(223);imagesc(Atau,Af,abs(SpecA')/max(max(abs(SpecA')))); axis([0 3 10 50]); axis xy
xlabel('Time (s)'); ylabel('Frequency (Hz)');

cbar = colorbar('horizontal'); % Create horizontal color bar
xlabel(cbar, 'Amplitude [•]'); % Label the color bar
caxis([0 1]);
set(cbar, 'Position', [0.36, 0.46, 0.2, 0.02]);

fidexmult = zeros(Num,length(Sig0)); clearvars Spec Af c bw IFfit extr_Sig findex fidexmult
num = Num;
for i = 1 : Num

    [Spec,Atau,Af] = mALCT(Sig0,Hz,Nh,Nc,num);
    Spec = Spec'; %
    c = findridges_m(Spec,deltaf,0.5,num,7);
    bw = 2;% 0.01*Hz;% the bandwidth of the TF filter for ridge extraction
    [IFfit,extr_Sig] = Dechirp_filter(Sig0,Hz,bw,Af(c),beta);
    fidexmult(i,:) = c;
    Sig0 = Sig0 - extr_Sig;
    num = num - 1;
end
% figure;plot(Af(fidexmult)');
% wfhowuifh
%% ridge path regrouping

Df = 0.05*Hz;
[findex,interset] = RPRG(fidexmult,Df);
% figure
% set(gcf,'Position',[20 100 640 500]);set(gcf,'Color','w');
% plot(t,Af(findex),'linewidth',1);
% xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
% ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
% set(gca,'FontSize',24);set(gca,'linewidth',2);axis([0 2 0 100]);

%% Smoothing of ridges

beta1 = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
iniIF = curvesmooth(Af(findex),beta1); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges

subplot(224);
plot(t,IF1,'k','LineWidth',2); hold on; plot(t,IF2,'b','LineWidth',2);
plot(t,iniIF(1,:),'--c','LineWidth',2);plot(t,iniIF(2,:),'--r','LineWidth',2); xlabel('Time (s)'); ylabel('Frequency (Hz)');axis([0 3 10 60]);
legend('T. IF1','T. IF2','E. IF1','E. IF2','Location','north','NumColumns',2)

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

% figure;,subplot(211);plot(t,iniIF(2,:)-IF1);subplot(212);plot(t,iniIF(1,:)-IF2);
%% run the CT with g for TFC domain
alpha1 = 50;
tt = -8:1/Hz:8;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0
eta = 1.5;

[tfc1, tfrtic, tcrtic] = sqSTCTm(Sig, eta/length(Sig), 1, h1);

% figure
% imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,201))), 1); axis xy; colormap(1-gray);
% xlabel('frequency (Hz)'); ylabel('chirp rate');
%% TFC ridges
if1 = iniIF(1,:);
if2 = iniIF(2,:);
chirpn1 = gradient(if1)/(t(2)-t(1));
chirpn2 = gradient(if2)/(t(2)-t(1));

line1 = zeros(length(Sig),Num); line2 = line1;
for i = 1:length(t)
    [~,bb] = min(abs(tfrtic-if1(i)/Hz));
    line1(i,1) = bb;
    [~,bb1] = min(abs(tcrtic-chirpn1(i)/Hz^2));
    line1(i,2) = bb1-1;
    [~,bb2] = min(abs(tfrtic-if2(i)/Hz));
    line2(i,1) = bb2;
    [~,bb3] = min(abs(tcrtic-chirpn2(i)/Hz^2));
    line2(i,2) = bb3-1;
end

% reconstruction
scrsz = get(0,'ScreenSize');
rrcon = zeros(2,length(Sig));
scale = alpha1;
for char = 1:length(Sig)
    tmp = zeros(2,2);
    tmp(1,1) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn1(char))));
    tmp(1,2) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn2(char))));
    tmp(2,1) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn1(char))));
    tmp(2,2) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn2(char))));

    xtmp = [tfc1(line1(char,2),line1(char,1),char); tfc1(line2(char,2),line2(char,1),char)]/Hz;
    rrcon(:,char) = tmp\xtmp;
end
%%
figure; set(gcf,'Position',[20 100 800 450]);
subplot(221)
plot(t,real(Sig1), 'color', 'k'); xlabel('time (s)');hold on
plot(t,real(rrcon(1,:)),'--c');ylim([-1.5 2]);
xlabel('Time (s)'); ylabel('Amplitude');legend('Real','Reconstructed','Orientation','horizontal');

subplot(223)
plot(t, real(Sig1')-real(rrcon(1,:)), 'r'); ylim([-1.5 2]); xlabel('Time (s)');ylabel('Amplitude');legend('Error')

subplot(222)
plot(t,real(Sig2), 'color', 'k'); xlabel('Time (s)');hold on;
plot(t,real(rrcon(2,:)),'--c');ylim([-1.5 2]);
xlabel('Time (s)'); ylabel('Amplitude');legend('Real','Reconstructed','Orientation','horizontal');

subplot(224)
plot(t, real(Sig2')-real(rrcon(2,:)), 'r'); ylim([-1.5 2]); xlabel('Time (s)');ylabel('Amplitude');legend('Error')

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
