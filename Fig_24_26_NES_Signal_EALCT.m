clear;clc;close all
load('AllData_Processed.mat')
% [(1:10)' Data.MaxForce']
Expi = 10;

dp = 25500;
x_o = Data.Disp(:,1,Expi); x = x_o(1:32:dp);             x = x'; 
t_o = Data.Time;  t = t_o(1:32:dp); Hz = 1/(t(2)-t(1));   t = t'; 

Sig = x'; % SampFreq = Hz;
figure; set(gcf,'Position',[20 100 800 450]);
subplot(221);plot(t_o,x_o,'k','linewidth',1); hold on; plot(t,x,'r--','linewidth',1);
legend('Original','Downsampled')
axis([0 1.5 -0.01 .01]); xlabel('Time (s)'); ylabel('Displacement (m)');

% -------------------- STFT ----------------------------
window = 128*4; Nfrebin = 1024;
subplot(222);
[Spec0,f] = STFT(real(Sig(:)),Hz,Nfrebin,window);
imagesc(t,f,abs(Spec0)/max(max(abs(Spec0)))); axis([0 1.5 0 50]);           

xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'YDir','normal');cbar = colorbar; ylabel(cbar, 'Amplitude [•]'); 
caxis([0 1]);
set(cbar, 'Position', [0.5, 0.2, 0.02, 0.33]); 

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

%% Estimating the number of components
Sigr = exp(2*pi*1i*(10*t));             %    IF =  10;
[Spec2,~] = STFT(real(Sigr(:)),Hz,Nfrebin,window);

N = length(t);
n = NaN(1,N);
hop = 8;                                                                   
deltatime = 23;
for i = deltatime :hop:N-(deltatime-1)
    Spec1 = 0*Spec0;
    Spec3 = 0*Spec2;
    pii = i-(deltatime-1); pf = i + (deltatime-1);
    Spec1(:,pii:pf) = Spec0(:,pii:pf);
    H1 = renyi(abs(Spec1));
    Spec3(:,pii:pf) = Spec2(:,pii:pf);
    H2 = renyi(abs(Spec3));
    n(i) = 2^(H1-H2);
end

figure;plot(floor(n),'bo')

%% ridge extraction and fitting
Num = max(floor(n));    % the number of the components
deltaf = floor(0.01*Hz);
beta = 1e-3;
Sig0 = real(Sig)';
if (isreal(Sig0))
    Sig0 = hilbert(Sig0);
end
fidexmult = zeros(Num,length(Sig0)); limits = zeros(Num,2); clearvars Spec 
num = Num;
% k = Num;
for i = 1 : Num

    Nh = 153; Nc = 2;
    [Spec,Atau,Af] = mALCT(Sig0,Hz,Nh,Nc,num);
    Spec = Spec'; 
    k = num;
    [c,lims] = findridges_m(Spec,deltaf,0.5,k,6);
    limits(i,:) = lims;
    bw = 17; %.028*Hz; % the bandwidth of the TF filter for ridge extraction
    [~,extr_Sig] = Dechirp_filter(Sig0,Hz,bw,Af(c),beta);
    fidexmult(i,:) = c;
    Sig0 = Sig0 - extr_Sig;
    num = num - 1;
end
figure;plot(t,Af(fidexmult)');


%% ridge path regrouping

Df = 0.01*Hz;
[findex,~] = RPRG(fidexmult,Df);
findex(findex<1)=1;
% figure
% set(gcf,'Position',[20 100 640 500]);set(gcf,'Color','w');
% plot(t,Af(findex),'linewidth',1);
% xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
% ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
% set(gca,'FontSize',24);set(gca,'linewidth',2);axis([0 1.5 0 45]);

%% Smoothing of ridges

beta1 = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
iniIF = curvesmooth(Af(findex),beta1); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges
% iniIF(iniIF<0)=1;

figure; set(gcf,'Position',[20 100 450 400]);
imagesc(t,f,abs(Spec0)); axis([0 1.5 0 50]);
hold on
plot(t(1:limits(1,2)),iniIF(1,1:limits(1,2)),'b','linewidth',2); hold on;
plot(t(1:limits(2,2)),iniIF(2,1:limits(2,2)),'r','linewidth',2); 
plot(t(1:limits(3,2)),iniIF(3,1:limits(3,2)),'color','[0, 0.5, 0]','linewidth',2);
% 
plot(t(limits(1,2)+1:end),iniIF(1,limits(1,2)+1:end),'--b','linewidth',2);
plot(t(limits(2,2)+1:end),iniIF(2,limits(2,2)+1:end),'--r','linewidth',2);
plot(t(limits(3,2)+1:end),iniIF(3,limits(3,2)+1:end),'color','[0, 0.5, 0]','LineStyle','--','linewidth',2);

% plot(t,iniIF);xlabel('Time / Sec'); ylabel('Frequency / Hz');
legend('NNM-1','NNM-2','Harmonic');axis([0 1.5 0 54]);

xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
axis xy
colormap(flipud(gray))

%% run the CT with g for TFC domain
alpha1 = 50;
tt = -8:1/Hz:8;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0

if (isreal(Sig))
    xm = hilbert(Sig);
end
eta = 2.0;
[tfc1, tfrtic, tcrtic] = sqSTCTm(xm, eta/length(xm), 1, h1);

% figure
% imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,201))), 1); axis xy; colormap(1-gray);
% xlabel('frequency (Hz)'); ylabel('chirp rate');
%% TFC ridges
if1 = iniIF(1,:); if2 = iniIF(2,:); if3 = iniIF(3,:); 
chirps = curvesmooth(gradient(iniIF)/(t(2)-t(1)),1e-4);
chirpn1 = chirps(1,:); chirpn2 = chirps(2,:); chirpn3 = chirps(3,:); 

line1 = zeros(length(Sig),2); line2 = line1; line3 = line1; 

for i = 1:length(xm)
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

%% reconstruction
% scrsz = get(0,'ScreenSize');
rrcon = zeros(Num,length(xm));
scale = alpha1;
stctmat = tfc1;
for char = 1:length(xm)
    tmp = zeros(Num,Num);
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
rrcon1 = rrcon; Sig1 = Sig;
%%
clc,clearvars -except rrcon1 Sig1 dp,%close all
load('AllData_Processed.mat')

Expi = 10;             

x_o = Data.Disp(:,2,Expi); x = x_o(1:32:dp);               x = x'; 
t_o = Data.Time;   t = t_o(1:32:dp); Hz = 1/(t(2)-t(1));   t = t'; 


Sig = x'; SampFreq = Hz;

figure; set(gcf,'Position',[20 100 800 450]);
subplot(221);plot(t_o,x_o,'k','linewidth',1); hold on; plot(t,x,'r--','linewidth',1);
axis([0 1.5 -0.01 .013]); xlabel('Time (s)'); ylabel('Displacement (m)');

%%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
window = 128*4; Nfrebin = 1024;
subplot(222);
[Spec0,f] = STFT(real(Sig(:)),Hz,Nfrebin,window);

imagesc(t,f,abs(Spec0)/max(max(abs(Spec0)))); axis([0 1.5 0 50]);           
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'YDir','normal');% set(gca,'FontSize',12);% set(gcf,'Color','w');
cbar = colorbar; % Create color bar
ylabel(cbar, 'Amplitude [•]'); % Label the color bar
caxis([0 1]);
set(cbar, 'Position', [0.45, 0.15, 0.02, 0.33]); % [left, bottom, width, height]                                                        % Added lims for Editor
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

%% Estimating the number of components
Sigr = exp(2*pi*1i*(10*t));               %  IF2 =  10;
[Spec2,~] = STFT(real(Sigr(:)),Hz,Nfrebin,window);

N = length(t);
n = NaN(1,N);
hop = 8;                                                                   
deltatime = 23;
for i = deltatime :hop:N-(deltatime-1)
    Spec1 = 0*Spec0;
    Spec3 = 0*Spec2;
    pii = i-(deltatime-1); pf = i + (deltatime-1);
    Spec1(:,pii:pf) = Spec0(:,pii:pf);
    H1 = renyi(abs(Spec1));
    Spec3(:,pii:pf) = Spec2(:,pii:pf);
    H2 = renyi(abs(Spec3));
    n(i) = 2^(H1-H2);
end

figure;plot(floor(n),'bo')
%% ridge extraction and fitting
Num = max(floor(n));    % the number of the components
deltaf = floor(0.01*Hz);
beta = 1e-3;
Sig0 = real(Sig)';
if (isreal(Sig0))
    Sig0 = hilbert(Sig0);
end
fidexmult = zeros(Num,length(Sig0)); limits = zeros(Num,2); clearvars Spec 
num = Num;
k = Num;
for i = 1 : Num

    Nh=153; Nc=2;
    [Spec,Atau,Af] = mALCT(Sig0,Hz,Nh,Nc,num);
    Spec = Spec'; %
    k = num;
    [c,lims] = findridges_m(Spec,deltaf,0.5,k,6);
    limits(i,:) = lims;
    bw = 17; %.028*Hz;% the bandwidth of the TF filter for ridge extraction
    [~,extr_Sig] = Dechirp_filter(Sig0,Hz,bw,Af(c),beta);
    fidexmult(i,:) = c;
    Sig0 = Sig0 - extr_Sig;
    num = num - 1;
end
% figure;plot(t,Af(fidexmult)');
% toughuitheg

%% ridge path regrouping

Df = 0.01*Hz;
[findex,interset] = RPRG(fidexmult,Df);
findex(findex<1)=1;
% figure
% set(gcf,'Position',[20 100 640 500]);set(gcf,'Color','w');
% plot(t,Af(findex),'linewidth',1);
% xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
% ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
% set(gca,'FontSize',24);set(gca,'linewidth',2);axis([0 1.5 0 45]);

%% Smoothing of ridges

beta1 = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
iniIF = curvesmooth(Af(findex),beta1); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges

figure; set(gcf,'Position',[20 100 450 400]);
imagesc(t,f,abs(Spec0)); axis([0 1.5 0 50]);
hold on
plot(t(1:limits(1,2)),iniIF(1,1:limits(1,2)),'b','linewidth',2); hold on;
plot(t(1:limits(2,2)),iniIF(2,1:limits(2,2)),'r','linewidth',2); 
plot(t(1:limits(3,2)),iniIF(3,1:limits(3,2)),'color','[0, 0.5, 0]','linewidth',2);
plot(t(1:limits(4,2)),iniIF(4,1:limits(4,2)),'color','[0.4940 0.1840 0.5560]','linewidth',2);
% 
plot(t(limits(1,2)+1:end),iniIF(1,limits(1,2)+1:end),'--b','linewidth',2);
plot(t(limits(2,2)+1:end),iniIF(2,limits(2,2)+1:end),'--r','linewidth',2);
plot(t(limits(3,2)+1:end),iniIF(3,limits(3,2)+1:end),'color','[0, 0.5, 0]','LineStyle','--','linewidth',2);
plot(t(limits(4,2)+1:end),iniIF(4,limits(4,2)+1:end),'color','[0.4940 0.1840 0.5560]','LineStyle','--','linewidth',2);

% plot(t,iniIF);xlabel('Time / Sec'); ylabel('Frequency / Hz');
legend('NNM-1','NNM-2','Harmonic','Harmonic');axis([0 1.5 0 54]);

xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
axis xy
colormap(flipud(gray))
% dgvg
%% run the CT with g for TFC domain
alpha1 = 50;
tt = -8:1/Hz:8;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0

if (isreal(Sig))
    xm = hilbert(Sig);
end
eta = 2.0;
[tfc1, tfrtic, tcrtic] = sqSTCTm(xm, eta/length(xm), 1, h1);

% figure
% imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,201))), 1); axis xy; colormap(1-gray);
% xlabel('frequency (Hz)'); ylabel('chirp rate');
%% TFC ridges
if1 = iniIF(1,:); if2 = iniIF(2,:); if3 = iniIF(3,:); if4 = iniIF(4,:);
chirps = curvesmooth(gradient(iniIF)/(t(2)-t(1)),1e-4);
chirpn1 = chirps(1,:); chirpn2 = chirps(2,:); chirpn3 = chirps(3,:); chirpn4 = chirps(4,:);

line1 = zeros(length(Sig),2); line2 = line1; line3 = line1; line4 = line1; 

for i = 1:length(xm)
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

    [~,bb6] = min(abs(tfrtic-if4(i)/Hz));
    line4(i,1) = bb6;
    [~,bb7] = min(abs(tcrtic-chirpn4(i)/Hz^2));
    line4(i,2) = bb7-1;
end

%% reconstruction
scrsz = get(0,'ScreenSize');
rrcon = zeros(Num,length(xm));
scale = alpha1;
stctmat = tfc1;
for char = 1:length(xm)
    tmp = zeros(Num,Num);
    tmp(1,1) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn1(char))));
    tmp(1,2) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn2(char))));
    tmp(1,3) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn3(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if3(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn3(char))));
    tmp(1,4) = 1/sqrt(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn4(char)))*exp(-pi*(tfrtic(line1(char,1))*Hz-if4(char))^2/(scale+1i*(tcrtic(line1(char,2))*Hz^2 - chirpn4(char))));

    tmp(2,1) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn1(char))));
    tmp(2,2) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn2(char))));
    tmp(2,3) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn3(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if3(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn3(char))));
    tmp(2,4) = 1/sqrt(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn4(char)))*exp(-pi*(tfrtic(line2(char,1))*Hz-if4(char))^2/(scale+1i*(tcrtic(line2(char,2))*Hz^2 - chirpn4(char))));

    tmp(3,1) = 1/sqrt(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line3(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn1(char))));
    tmp(3,2) = 1/sqrt(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line3(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn2(char))));
    tmp(3,3) = 1/sqrt(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn3(char)))*exp(-pi*(tfrtic(line3(char,1))*Hz-if3(char))^2/(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn3(char))));
    tmp(3,4) = 1/sqrt(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn4(char)))*exp(-pi*(tfrtic(line3(char,1))*Hz-if4(char))^2/(scale+1i*(tcrtic(line3(char,2))*Hz^2 - chirpn4(char))));

    tmp(4,1) = 1/sqrt(scale+1i*(tcrtic(line4(char,2))*Hz^2 - chirpn1(char)))*exp(-pi*(tfrtic(line4(char,1))*Hz-if1(char))^2/(scale+1i*(tcrtic(line4(char,2))*Hz^2 - chirpn1(char))));
    tmp(4,2) = 1/sqrt(scale+1i*(tcrtic(line4(char,2))*Hz^2 - chirpn2(char)))*exp(-pi*(tfrtic(line4(char,1))*Hz-if2(char))^2/(scale+1i*(tcrtic(line4(char,2))*Hz^2 - chirpn2(char))));
    tmp(4,3) = 1/sqrt(scale+1i*(tcrtic(line4(char,2))*Hz^2 - chirpn3(char)))*exp(-pi*(tfrtic(line4(char,1))*Hz-if3(char))^2/(scale+1i*(tcrtic(line4(char,2))*Hz^2 - chirpn3(char))));
    tmp(4,4) = 1/sqrt(scale+1i*(tcrtic(line4(char,2))*Hz^2 - chirpn4(char)))*exp(-pi*(tfrtic(line4(char,1))*Hz-if4(char))^2/(scale+1i*(tcrtic(line4(char,2))*Hz^2 - chirpn4(char))));

    xtmp = [stctmat(line1(char,2),line1(char,1),char); stctmat(line2(char,2),line2(char,1),char); stctmat(line3(char,2),line3(char,1),char); stctmat(line4(char,2),line4(char,1),char)]/Hz;
    rrcon(:,char) = tmp\xtmp;
end
%%
figure;
rrcon1(1,1:limits1(1,1)) = zeros(1,limits1(1,1)); rrcon1(1,limits1(1,2):length(Sig)) = zeros(1,length(Sig)-limits1(1,2)+1);
rrcon(1,1:limits(1,1)) = zeros(1,limits(1,1)); rrcon(1,limits(1,2):length(Sig)) = zeros(1,length(Sig)-limits(1,2)+1);

subplot(422);plot(t,real(rrcon(1,:)),'k',t,real(rrcon1(1,:)),'b'); ylabel('Displacement (m)'); 
title('NNM-1'); legend('NES','LO'); axis([0 1.5 -0.008 0.008]);

rrcon1(2,1:limits1(2,1)) = zeros(1,limits1(2,1)); rrcon1(2,limits1(2,2):length(Sig)) = zeros(1,length(Sig)-limits1(2,2)+1);
rrcon(2,1:limits(2,1)) = zeros(1,limits(2,1)); rrcon(2,limits(2,2):length(Sig)) = zeros(1,length(Sig)-limits(2,2)+1);

subplot(423);plot(t,real(rrcon(2,:)),'k',t,real(rrcon1(2,:)),'b');ylabel('Displacement (m)'); title('NNM-2');set(gca,'fontsize',10);axis([0 1.5 -0.004 0.004]);
title('NNM-2'); legend('NES','LO'); axis([0 1.5 -0.004 0.004]);

rrcon1(3,1:limits1(3,1)) = zeros(1,limits1(3,1)); rrcon1(3,limits1(3,2):length(Sig)) = zeros(1,length(Sig)-limits1(3,2)+1);
rrcon(3,1:limits(3,1)) = zeros(1,limits(3,1)); rrcon(3,limits(3,2):length(Sig)) = zeros(1,length(Sig)-limits(3,2)+1);

subplot(421);plot(t,real(rrcon(3,:)),'k',t,real(rrcon1(3,:)),'b');ylabel('Displacement (m)'); title('Harmonic'); set(gca,'fontsize',10);axis([0 1.5 -0.002 0.002]);
title('Harmonic'); legend('NES','LO'); axis([0 1.5 -0.003 0.003]);

rrcon(4,1:limits(4,1)) = zeros(1,limits(4,1)); rrcon(4,limits(4,2):length(Sig)) = zeros(1,length(Sig)-limits(4,2)+1);
subplot(424);plot(t,real(rrcon(4,:)),'k');ylabel('Displacement (m)'); title('Harmonic');  set(gca,'fontsize',10);axis([0 1.5 -0.0011 0.0011]);
title('Harmonic'); legend('NES'); axis([0 1.5 -0.0011 0.0011]);

subplot(425);plot(t,Sig,'k');hold on;plot(t,real(sum(rrcon)),'--c');title('NES');legend('Real','Reconstructed');
ylabel('Displacement (m)');axis([0 1.5 -0.01 0.013]);
subplot(426);plot(t,real(sum(rrcon))-Sig','r');legend('Error');axis([0 1.5 -0.01 0.013]);ylabel('Displacement (m)');

subplot(427);plot(t,Sig1,'k');hold on;plot(t,real(sum(rrcon1)),'--c');title('LO');legend('Real','Reconstructed');
xlabel('Time (s)');  ylabel('Displacement (m)');axis([0 1.5 -0.01 0.013]);
subplot(428);plot(t,real(sum(rrcon1))-Sig1','r');legend('Error');axis([0 1.5 -0.01 0.013]);xlabel('Time (s)');  ylabel('Displacement (m)');

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
