clc; clear;  close all
Hz = 100;
t = 0:1/Hz:3; t = t';

am1 = ones(1,length(t)); am2 = 0.9*am1;

figure; set(gcf,'Position',[20 100 800 450]);
subplot(221);plot(t,am2,'color','0.25 0.80 0.54','LineWidth',2);hold on;plot(t,am1,'color','0.5 0 1','LineWidth',2);
axis xy; xlim([0,t(end)]);ylim([0 2]);
legend('Mode 1','Mode 2','Location','northwest')
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');

Sig1 = exp(2*pi*1i*(30*t + 3*sin(3*t)));            IF2 = 30 + 9*cos(3*t);
Sig2 = 0.9*exp(2*pi*1i*(30*t - 3*sin(3*t) + 3/(2*pi)));            IF1 = 30 - 9*cos(3*t);

initstate(2) ;
x = Sig1 + Sig2; x = x';
x=awgn(x,10);

%%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
window = 128/2; Nfrebin = 1024;
figure;
[Spec0,f] = STFT(real(x(:)),Hz,Nfrebin,window);
imagesc(t,f,abs(Spec0)); axis([0 2 0 100]);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'YDir','normal');% set(gca,'FontSize',12);% set(gcf,'Color','w');

% x = awgn(x,10);
threshold = 0.03;

ng = 151;
Nx = length(x);
nw_2 = Nx;
ntau = Nx;
nc = 71;
% figure; plot(t,f1);hold on;plot(t,f2);plot(t,f3);
% xlim([0,tend]); ylim([0,fs/2]);
% xlabel('Time (s)'); ylabel('Freq (Hz)');
% [TFR,Atau,Af] = ALCT(x,fs,Nh,Ntau,Nfext,Nc,threshold);

% figure;imagesc(Atau,Af,abs(TFR.'));axis xy;
% xlim([0,tend]);ylim([0,fs/2]);
% xlabel('Time (s)');ylabel('Freq (Hz)');
% [TFR2] = TFRnremove(TFR,0.05);
% figure;imagesc(Atau,Af,abs(TFR2.'));axis xy;
% xlim([0,tend]);ylim([0,fs/2]);
% xlabel('Time (s)');ylabel('Freq (Hz)');

% function [TFR,Atau,Af] = ALCT(x,fs,Nh,Ntau,Nfext,c_n,threshold)
% L : length of discrete signal
% fs : sampling frequency
% ntau : LCT time
% nw : linear chirplet transform (LCT) frequency
% nc : LCT chirp rate
% ng : LCT chirp window.
% nw_2 and nc_2 are about two times of nw and nc.

L = length(x);
nw = round(nw_2/2);
nc = nc+1-rem(nc,2);
ng = ng+1-rem(ng,2);
ng_h = (ng-1)/2;
nc_2 = nc*2-1;
cend = Hz*Hz/ng;

% Atau, Aw and Ac are arguments of LCT time, frequency and chirp rate, respectively.
Atau = (1/ntau:1/ntau:1)*L/Hz;
Anw_2 = (0:(nw_2-1))*Hz/nw_2;
Aw = Anw_2(1:nw);
Ac = (-1:2/(nc-1):1)*cend;
% Ac = -65:1:65;                          %%%%%%% modified by CL  %%%%%%%
%%
% Calculation of signal LCT
g = 0.54 - 0.46*cos(2.0*pi*(1:ng)'/(ng+1));  % HAMMING window
g = g/sum(g);
taon = round(Atau*Hz);
x = hilbert(x);
xext = [zeros(1,ng_h) x zeros(1,ng_h) ];     % extended signal
xmat = xext(repmat((0:ng-1),[ntau,1,nc])+repmat(taon.',[1,ng,nc]));% matrix of truncated signal
gmat = repmat(g.',[ntau,1,nc]);
gt = (-ng_h:ng_h)/Hz;                        % time of window
kn(1,:,:) = exp(-1j*pi*(gt).^2.'*Ac);        % kernel of LCT
knmat = repmat(kn,[ntau,1,1]);
xLCTa = xmat.*gmat.*knmat;
xLCTb = zeros(ntau,nw_2,nc);
xLCTb(:,1:ng,:) = xLCTa;
xLCTc = circshift(xLCTb,-ng_h,2);
xLCTc = fft(xLCTc,[],2);
xLCT = xLCTc(:,1:nw,:);

% Calculation of LCT Reference
Ac_2 = (-2:2/(nc-1):2)*cend;
% Ac_2 = -130:1:130;                          %%%%%%% modified by CL  %%%%%%%

RefLCTb = zeros(nw_2,nc_2);
Refhmat = repmat(g,[1,nc_2]);
Refknmat = exp(-1j*pi*(gt).^2.'*Ac_2);
RefLCTa = Refhmat.*Refknmat;
RefLCTb(1:ng,:) = RefLCTa;
RefLCTc = circshift(RefLCTb,-ng_h,1);
RefLCTc = fft(RefLCTc,[],1);
LCTRef = fftshift(RefLCTc,1);

fmidn = round((nw_2+1)/2)*ones(1,ntau); % frequency center index of LCT Reference
cmidn = (nc_2+1)/2*ones(1,ntau);        % chirp rate center index of LCT Reference

%iterative process
TFR = zeros(ntau,nw);                       %% ---- for the other component, run this line again ---- %%% 
mf = 0 : nw-1;
mc(1,1,:) = (0:nc-1)*nw_2;
mfc = repmat(mf,ntau,1,nc)+repmat(mc,ntau,nw,1);
% xenergy = sum(sum(xLCT(:,:,(nc+1)/2).*conj(xLCT(:,:,(nc+1)/2))));%signal energy
xenergy = sum(sum(xLCT(:,:,(nc+0)/1).*conj(xLCT(:,:,(nc+0)/1))));%signal energy

% ratio = 1;
%% while ratio>threshold
% Step 2
absxLCT = abs(xLCT);
[~,TFCMaxIndex] = max(absxLCT,[],[2,3],'linear');
[~,maxfn,maxcn] = ind2sub([ntau,nw,nc],TFCMaxIndex);   % positions in matrix % CL (3,1);, (2,2);, ...
TFMaxIndex = sub2ind([ntau,nw],1:ntau,maxfn.');
mode = xLCT(TFCMaxIndex.');
TFR(TFMaxIndex) = TFR(TFMaxIndex) + mode;

% Step 3
fdist = fmidn.'- maxfn + 1;                  % distance of frequency
cdist = cmidn.'- maxcn + 1;                  % distance of chirp rate
mtau = fdist + (cdist-1)*nw_2;
m1b = repmat(mtau,1,nw,nc) + mfc;
xLCTk = repmat(mode.',1,nw,nc).*LCTRef(m1b);  % LCT of the kth mode
% Step 4
xLCT = xLCT-xLCTk;% residual LCT
% Step 5
renergy = sum(sum(xLCT(:,:,(nc+1)/2).*conj(xLCT(:,:,(nc+1)/2))));%residual energy
ratio = renergy/xenergy;
% end
% end
figure; set(gcf,'Position',[20 100 900 550]);
subplot(232);imagesc(Atau,Aw,abs(TFR.'));axis xy;
xlim([0,t(end)]);ylim([0 100])
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',12,'FontName','Times New Roman');
colormap(flipud(gray))
% ax.CLim([0 1.2])

subplot(231);plot(t,IF1,'b');hold on;plot(t,IF2,'r');axis xy;
xlim([0,t(end)]);ylim([0 100]);legend('Mode 1','Mode 2','Location','north')
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',12,'FontName','Times New Roman');

TFR = zeros(ntau,nw);   
absxLCT = abs(xLCT);
[~,TFCMaxIndex] = max(absxLCT,[],[2,3],'linear');
[~,maxfn,maxcn] = ind2sub([ntau,nw,nc],TFCMaxIndex);   % positions in matrix % CL (3,1);, (2,2);, ...
TFMaxIndex = sub2ind([ntau,nw],1:ntau,maxfn.');
mode2 = xLCT(TFCMaxIndex.');
TFR(TFMaxIndex) = TFR(TFMaxIndex) + mode2;
% Step 3
fdist = fmidn.'- maxfn + 1;                  % distance of frequency
cdist = cmidn.'- maxcn + 1;                  % distance of chirp rate
mtau = fdist + (cdist-1)*nw_2;
m1b = repmat(mtau,1,nw,nc) + mfc;
xLCTk = repmat(mode2.',1,nw,nc).*LCTRef(m1b);  % LCT of the kth mode
% Step 4
xLCT = xLCT-xLCTk;% residual LCT
% Step 5
renergy = sum(sum(xLCT(:,:,(nc+1)/2).*conj(xLCT(:,:,(nc+1)/2))));%residual energy
ratio = renergy/xenergy;

subplot(233);imagesc(Atau,Aw,abs(TFR.'));axis xy;
xlim([0,t(end)]);ylim([0 100])
xlabel('Time (s)','FontSize',12,'FontName','Times New Roman');
ylabel('Frequency (Hz)','FontSize',12,'FontName','Times New Roman');
colormap(flipud(gray))
% clim([0 1.2])
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

figure; 
set(gcf,'Position',[20 100 800 450]);	
subplot(221);plot(t,real(Sig1),'k');hold on;plot(t,real(mode),'--c');ylim([-1.5 2])
xlabel('Time (s)');ylabel('Amplitude');legend('Original','Reconstructed','Orientation','horizontal');
subplot(223);plot(t,real(Sig1'-mode),'r');ylim([-1.5 2]);xlabel('Time (s)');ylabel('Amplitude');legend('Error')
subplot(222);plot(t,real(Sig2),'k');hold on;plot(t,real(mode2),'--c');ylim([-1.5 2])
xlabel('Time (s)');ylabel('Amplitude');legend('Original','Reconstructed','Orientation','horizontal');
subplot(224);plot(t,real(Sig2'-mode2),'r');xlabel('Time (s)');ylabel('Amplitude');
ylim([-1.5 2]);xlabel('Time (s)');ylabel('Amplitude');legend('Error')

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
