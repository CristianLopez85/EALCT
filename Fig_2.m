clc; clear; close all
Hz = 400;
t = 0:1/Hz:2; t = t';

Sig1 = exp(2*pi*1i*(-20*t.^2 + 90*t));               IF1 = -40*t + 90;
Sig2 = exp(2*pi*1i*(20*t.^2 + 10*t));                IF2 =  40*t + 10;
Sig3 = exp(2*pi*1i*(30*t));                          IF3 = 30*ones(length(t),1);

Sig1 = [zeros(100,1); Sig1(101:700); zeros(101,1)];
Sig3 = [Sig3(1:400); zeros(401,1)];

Sig = Sig1 + Sig2 + Sig3;

%% Plot time series
figure; set(gcf,'Position',[200 100 720 200]);
subplot(133)
plot(t,real(Sig1), 'color', 'k'); ylim([-1.5 1.5]); 
xlabel('Time (s)'); ylabel('Amplitude');

subplot(131)
plot(t,real(Sig2), 'color', 'k');ylim([-1.5 1.5]); 
xlabel('Time (s)'); ylabel('Amplitude');

subplot(132)
plot(t,real(Sig3), 'color', 'k'); ylim([-1.5 1.5]); 
xlabel('Time (s)'); ylabel('Amplitude');
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
%%

%%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
window = 128; Nfrebin = 1024;
% figure
[Spec0,f] = STFT(real(Sig(:)),Hz,Nfrebin,window);
% imagesc(t,f,abs(Spec0)); axis([0 2 0 100]);
% set(gcf,'Position',[20 100 320 250]);
% xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
% ylabel('Frequency / Hz','FontSize',12,'FontName','Times New Roman');
% set(gca,'YDir','normal');

% --------- Estimating the number of components ---------------
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
% -------------- Fig 3 --------------------------- %%
subplot(222);plot(t,floor(n),'b.')
xlabel('Time (s)'); ylabel('Number of components'); ylim([0 3])
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

Num = max(floor(n));    % the number of the components

%% ridge extraction and fitting
delta = 8;
beta = 1e-3;
Sig0 = real(Sig)';
if (isreal(Sig0))
    Sig0 = hilbert(Sig0);
end
fidexmult = zeros(Num,length(Sig0)); limits = zeros(Num,2); clearvars Spec Af c bw IFfit extr_Sig findex fidexmult
num = Num;
Nh=251; Nc=71; % threshold=0.03; 
for i = 1 : Num
    
    [Spec,Atau,Af] = mALCT(Sig0,Hz,Nh,Nc,num);
    Spec = Spec'; 

    [c,lims] = findridges_m(Spec,delta,0.5,num,20);
    limits(i,:) = lims;
    c(c<1)=num;
    bw = t(end)/100;    % the bandwidth of the TF filter for ridge extraction
    [IFfit,extr_Sig] = Dechirp_filter(Sig0,Hz,bw,Af(c),beta);
    fidexmult(i,:) = c;
    Sig0 = Sig0 - extr_Sig;
    num = num - 1;
end

figure;plot(Af(fidexmult)')

%% ridge path regrouping

Df = 0.03*Hz;
[findex,interset] = RPRG(fidexmult,Df);
figure
set(gcf,'Position',[20 100 640 500]);set(gcf,'Color','w');
plot(t,Af(findex),'linewidth',1);
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency / Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);set(gca,'linewidth',2);axis([0 2 0 100]);

% Original and new ridges
org_ridg = Af(fidexmult); 
new_ridg = Af(findex);

% Ridges comparison
result = org_ridg==new_ridg;

for i = 1 : Num
    for j = 1 : 2
        if 1 ~= result(i,limits(i,j))
            lim_ch(i,j) = limits(i,j);
        end
    end 
end
[aa,~] = find(lim_ch); aa1 = sort(aa,'descend');
limits(aa) = limits(aa1); 

%% Smoothing of ridges

beta = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
iniIF = curvesmooth(Af(findex),beta); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges
% 
figure
plot(t,iniIF);xlabel('Time / Sec'); ylabel('Frequency / Hz');axis([0 2 0 100]);
hold on
plot(t,IF1,'linewidth',1); plot(t,IF2,'linewidth',1); plot(t,IF3,'linewidth',1);
% ucvkc
% figure;,subplot(211);plot(t,iniIF(2,:)-IF1);subplot(212);plot(t,iniIF(1,:)-IF2);

% Fig 1a
figure; set(gcf,'Position',[20 100 320 250]);
plot(t(101:700),IF1(101:700),'Color','b','linewidth',2); hold on;
plot(t,IF2,'color','r','linewidth',2);
plot(t(1:400),IF3(1:400),'color','[0, 0.5, 0]','linewidth',2); 
xlabel('Time (s)'); ylabel('Frequency (Hz)');
legend('Mode 1','Mode 2','Mode 3')
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

% Fig 1b
figure; set(gcf,'Position',[20 100 320 250]);
plot(t(limits(3,1):limits(3,2)),iniIF(3,limits(3,1):limits(3,2)),'color','r','linewidth',2); hold on;
plot(t(limits(2,1):limits(2,2)),iniIF(2,limits(2,1):limits(2,2)),'color','b','linewidth',2); 
plot(t(limits(1,1):limits(1,2)),iniIF(1,limits(1,1):limits(1,2)),'Color','[0, 0.5, 0]','linewidth',2); 

plot(t(1:limits(2,1)-1),iniIF(2,1:limits(2,1)-1),'--b','linestyle','--','LineWidth',1); 
plot(t(limits(2,1):100),iniIF(2,limits(2,1):100),'--b','linestyle','--','LineWidth',1);  
plot(t(limits(2,2)+1:end),iniIF(2,limits(2,2)+1:end),'--b','linestyle','--','LineWidth',1);
plot(t(limits(1,2)+1:end),iniIF(1,limits(1,2)+1:end),'Color','[0, 0.5, 0]','linestyle','--','LineWidth',1);

xlabel('Time (s)'); ylabel('Frequency (Hz)');
% legend('Real IF1','Real IF2','Real IF3','Location','north')
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

%% run the CT with g for TFC domain
alpha1 = 50;
tt = -8:1/Hz:8;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0

xm = Sig;
[tfc1, tfrtic, tcrtic] = sqSTCTm(xm, 1.5/length(xm), 1, h1);

figure
imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,201))), 1); axis xy; colormap(1-gray);
xlabel('frequency (Hz)'); ylabel('chirp rate');

%% TFC ridges

if1 = iniIF(1,:);
if2 = iniIF(2,:);
if3 = iniIF(3,:);
chirpn1 = gradient(if1)/(t(2)-t(1));
chirpn2 = gradient(if2)/(t(2)-t(1));
chirpn3 = gradient(if3)/(t(2)-t(1));

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

% reconstruction
scrsz = get(0,'ScreenSize');
rrcon = zeros(3,length(xm));
scale = alpha1;
stctmat = tfc1;
for char = 1:length(xm)
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

figure;
subplot(233)
plot(t,real(Sig1), 'color', 'k'); xlabel('time (s)');hold on
rrcon(2,1:limits(2,1)) = zeros(1,limits(2,1)); rrcon(2,limits(2,2):length(Sig)) = zeros(1,length(Sig)-limits(2,2)+1);
plot(t,real(rrcon(2,:)),'--c');ylim([-1.5 1.5]); 
xlabel('Time (s)'); ylabel('Amplitude');%legend('Real','Reconstructed','Orientation','horizontal'); set(gca,'fontsize',10)

subplot(236)
plot(t, real(Sig1')-real(rrcon(2,:)), 'r'); ylim([-1.5 1.5]);xlabel('Time (s)');ylabel('Amplitude');set(gca,'fontsize',10);%legend('Error')

subplot(231)
plot(t,real(Sig2), 'color', 'k'); xlabel('Time (s)');hold on;
plot(t,real(rrcon(3,:)),'--c');ylim([-1.5 1.5]); 
xlabel('Time (s)'); ylabel('Amplitude');%legend('Real','Reconstructed','Orientation','horizontal'); set(gca,'fontsize',10)

subplot(234)
plot(t, real(Sig2')-real(rrcon(3,:)), 'r'); ylim([-1.5 1.5]);xlabel('Time (s)');ylabel('Amplitude');set(gca,'fontsize',10);%legend('Error')

subplot(232)
plot(t,real(Sig3), 'color', 'k'); xlabel('Time (s)');hold on;
rrcon(1,limits(1,2):length(Sig)) = zeros(1,length(Sig)-limits(1,2)+1);
plot(t,real(rrcon(1,:)),'--c');ylim([-1.5 1.5]); 
xlabel('Time (s)'); ylabel('Amplitude');legend('Real','Reconstructed','Orientation','horizontal'); set(gca,'fontsize',10)

subplot(235)
plot(t, real(Sig3')-real(rrcon(1,:)), 'r'); ylim([-1.5 1.5]);xlabel('Time (s)');ylabel('Amplitude');set(gca,'fontsize',10);legend('Error')
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')