%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Copyright 2016 Boualem Boashash
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Authors: Prof. Boualem Boashash   (boualem.boashash@gmail.com)
%          Samir Ouelha             (samir_ouelha@hotmail.fr)
%          
% The following references should be cited whenever this script is used:
% [1] B. Boashash, Samir Ouelha, Designing high-resolution time-frequency
% and time-scale distributions for the analysis and  classification of non 
% stationary signals: a tutorial review with a comparison of features performance
%, Digital Signal Processing, In Press, 2017.
% [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package to compute
% Time-Frequency Distributions and related Time-Scale methods with extraction
% of signal characteristics, SoftwareX, In Press, 2017.
%
% This study was funded by grants from the ARC and QNRF NPRP 6-885-2-364
% and NPRP 4-1303-2-517 
%
%
%%
clear;clc;close all
% load EEG_examples
load seizure
Hz = 16; t = 0:1/Hz:8-1/Hz;
Sig = Sig';
% figure
% set(gcf,'Position',[20 100 320 250]); set(gcf,'Color','w');
% plot(t,Sig,'linewidth',1);
% xlabel('Time / Sec','FontSize',12,'FontName','Times New Roman');
% ylabel('Amplitude','FontSize',12,'FontName','Times New Roman');
% set(gca,'FontSize',12);% set(gca,'linewidth',1);

%%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
window = 128/2; Nfrebin = 1024;
[Spec0,f] = STFT(real(Sig(:)),Hz,Nfrebin,window);
figure; set(gcf,'Position',[20 100 800 250]);
subplot(121);imagesc(t,f,abs(Spec0)); axis([0 t(end) 0 Hz/2]);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
set(gca,'YDir','normal'); set(gcf,'Color','w');

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

%  figure;plot(floor(n),'bo')
%% ridge extraction and fitting
Num = max(floor(n));    % the number of the components
delta = 20; %0.1875*16;
beta = 1e-3;
Sig0 = real(Sig)';
if (isreal(Sig0))
    Sig0 = hilbert(Sig0);
end
fidexmult = zeros(Num,length(Sig0)); limits = zeros(Num,2); clearvars Spec
num = Num;
for i = 1 : Num

    Nh = 128; Nc= 80; % threshold=0.03;
    [Spec,Atau,Af] = mALCT(Sig0,Hz,Nh,Nc,num);
    Spec = Spec'; 

    [c,lims] = findridges_m(Spec,delta,0.5,num,7);
    bw = Hz/60; % .01*Hz;% the bandwidth of the TF filter for ridge extraction
    [~,extr_Sig] = Dechirp_filter(Sig0,Hz,bw,Af(c),beta);
    fidexmult(i,:) = c;
    Sig0 = Sig0 - extr_Sig;
    num = num - 1;
end

%% ridge path regrouping

Df = length(Af)/15; %.03*Hz;
[findex,interset] = RPRG(fidexmult,Df);

%% Smoothing of ridges

beta1 = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
iniIF = curvesmooth(Af(findex),beta1); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges

subplot(122)
plot(t,iniIF(1,:),'b','linewidth',2); hold on;
plot(t,iniIF(2,:),'r','linewidth',2);plot(t,iniIF(3,:),'color','[0, 0.5, 0]','linewidth',2);
xlabel('Time / Sec'); ylabel('Frequency / Hz');axis([0 t(end) 0 8]);
legend('Mode 1','Mode 2','Mode 3','Mode 4','Location','Northeast')

set(findall(gcf,'-property','FontSize'),'FontSize',14, 'FontName', 'Times New Roman')

%% run the CT with g for TFC domain
alpha1 = 50;
tt = -8:1/Hz:8;
h1 = exp(-pi*alpha1*tt.^2'); % window g_0

if (isreal(Sig))
    xm = hilbert(Sig);
end
eta = 1;
[tfc1, tfrtic, tcrtic] = sqSTCTm(xm, eta/length(xm), 1, h1);

% figure
% imageSQ(Hz*tfrtic, Hz^2*tcrtic, abs(squeeze(tfc1(:,:,201))), 1); axis xy; colormap(1-gray);
% xlabel('frequency (Hz)'); ylabel('chirp rate');
%% TFC ridges
if1 = iniIF(1,:); if2 = iniIF(2,:); if3 = iniIF(3,:); % if4 = iniIF(4,:);
chirpn1 = gradient(if1)/(t(2)-t(1)); chirpn1 = mean(chirpn1)*ones(1,length(Sig));
chirpn2 = gradient(if2)/(t(2)-t(1)); chirpn2 = mean(chirpn2)*ones(1,length(Sig));
chirpn3 = gradient(if3)/(t(2)-t(1)); chirpn3 = mean(chirpn3)*ones(1,length(Sig));

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
figure; set(gcf,'Position',[20 100 800 450]);
subplot(321);plot(t,real(rrcon(1,:)),'k');xlabel('Time (s)'); ylabel('Amplitude'); set(gca,'fontsize',10); axis([0 t(end) -0.12 0.12]);
subplot(322);plot(t,real(rrcon(2,:)),'k');xlabel('Time (s)'); ylabel('Amplitude'); set(gca,'fontsize',10); axis([0 8 -0.1 0.1]);
subplot(323);plot(t,real(rrcon(3,:)),'k');xlabel('Time (s)'); ylabel('Amplitude'); set(gca,'fontsize',10); axis([0 8 -0.07 0.07]);

subplot(325);plot(t,Sig,'k');hold on;plot(t,real(sum(rrcon)),'--c');axis([0 8 -0.15 0.2]); xlabel('Time (s)'); ylabel('Amplitude'); 
legend('Real','Reconstructed','Orientation','horizontal');
subplot(326);plot(t,real(sum(rrcon))-Sig','r');axis([0 8 -0.15 0.2]);xlabel('Time (s)'); ylabel('Amplitude'); legend('Error');

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')