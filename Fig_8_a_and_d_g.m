clear;clc;close all
% profile on
tic
warning('off')
SampFreq = 400;
t = 0:1/SampFreq:2; t = t';

initstate(1) ;
am1 = smooth(cumsum(randn(length(t),1)) ./ SampFreq, 200, 'loess') ;
am1 = 2 + am1 ./ max(abs(am1)) ;
am2 = smooth(cumsum(randn(length(t),1)) ./ SampFreq, 200, 'loess') ;
am2 = 2 + am2 ./ max(abs(am2)) ;

Sig1 = am1.*exp(2*pi*1i*(-20*t.^2 + 90*t));               IF1 = -40*t + 90;
Sig2 = am2.*exp(2*pi*1i*(20*t.^2 + 10*t));                IF2 =  40*t + 10;

Sig = Sig1 + Sig2;
Sig = Sig';

figure; set(gcf,'Position',[20 100 800 450]); subplot(221)
plot(t,Sig,'k');xlabel('Time (s)'); ylabel('Amplitude');

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

window = 128; Nfrebin = 1024;

[m,n] = size(Sig);
time = (1:n)/SampFreq;
fre = (SampFreq/2)/(n/2):(SampFreq/2)/(n/2):(SampFreq/2);

x1 = 0.7;  x2 = 1.3;
y1 = 30;   y2 = 70;

lw = 200;

%......................... Get SET and SST ...........................
[Ts] = SST_Y(Sig',lw);
[~, Te, ~] = SET_Y(Sig',lw);

figure; set(gcf,'Position',[20 100 800 450]);

subplot(221);imagesc(time,fre,abs(Ts)/max(max(abs(Ts)))); ylim([0 110])
xlabel('Time (s)');ylabel('Frequency (Hz)');axis xy

subplot(222);imagesc(time,fre,abs(Te)/max(max(abs(Te))));
xlabel('Time (s)');ylabel('Frequency (Hz)');axis xy; ylim([0 110])

%......................... Get GLCT ...........................

N = 7;
hlength = 100;
tfr = GLCT((real(Sig))',N,SampFreq,hlength);
subplot(223);imagesc(time,fre,abs(tfr).^2/max(max(abs(tfr).^2)));
axis xy
ylabel('Frequency (Hz)'); xlabel('Time (s)');  ylim([0 110])

%......................... Get ACMD ...........................
Sig = real(Sig);
% Fourier spectrum
Spec = 2*abs(fft(Sig))/length(Sig);
Spec = Spec(1:round(end/2));
Freqbin = linspace(0,SampFreq/2,length(Spec));
	
% parameter setting
alpha0 = 1e-3; % or 1e-3  %or alpha = 2e-4; % if this parameter is larger, it will help the algorithm to find correct modes even the initial IFs are too rough. But it will introduce more noise and also may increase the interference between the signal modes
beta = 1e-4; % 0.5e-4  this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
tol = 1e-8;%
% Component 1 extraction
[~,findex1] = max(Spec);
peakfre1 = Freqbin(findex1); % IF initialization by finding peak frequency of the Fourier spectrum
iniIF1 = peakfre1*ones(1,length(Sig)); % initial IF vector

[IFest1 IAest1 sest1] = ACMD(Sig,SampFreq,iniIF1,alpha0,beta,tol); % extract the signal component and its IA and IF
% Component 2 extraction
residue1 = Sig - sest1; % obtain the residual signal by removing the extracted component from the raw signal
resSpec = 2*abs(fft(residue1))/length(residue1); % FFT of the residual signal
resSpec = resSpec(1:round(end/2));
[~,findex2] = max(resSpec);
peakfre2 = Freqbin(findex2); % IF initialization by finding peak frequency of the Fourier spectrum of the residue signal
iniIF2 = peakfre2*ones(1,length(Sig)); % initial IF vector

[IFest2 IAest2 sest2] = ACMD(residue1,SampFreq,iniIF2,alpha0,beta,tol); % extract the component 2 and its IA and IF from the residue signal

% adaptive time-frequency spectrum
band = [0 SampFreq/2];
IFmulti = [IFest1;IFest2];
IAmulti = [IAest1;IAest2];
[ASpec fbin] = TFspec(IFmulti,IAmulti,band);
subplot(224);imagesc(t,fbin,abs(ASpec)/max(max(abs(ASpec)))); 	 
xlabel('Time (s)');ylabel('Frequency (Hz)'); ylim([0 110])
set(gca,'YDir','normal')	


set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')

function [IF, Te, tfr] = SET_Y(x,hlength)
%   Synchroextracting Transform
%	x       : Signal.
%	hlength : Window length.

%   IF   : Synchroextracting operator representation.
%   Te   : SET result.
%   tfr  : STFT result
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%
%   Written by YuGang in Shandong University at 2016.5.13.

[xrow,xcol] = size(x);

N=xrow;

if (xcol~=1)
    error('X must be column vector');
end

if (nargin < 2)
    hlength=round(xrow/8);
end

t=1:N;
% ft = 1:round(N/2);

[~,tcol] = size(t);

hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);ht=ht';

% Gaussian window
h = exp(-pi/0.32^2*ht.^2);
% derivative of window
dh = -2*pi/0.32^2*ht .* h; % g'

[hrow,~]=size(h); Lh=(hrow-1)/2;

tfr1= zeros (N,tcol);
tfr2= zeros (N,tcol);

for icol=1:tcol
    ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices= rem(N+tau,N)+1;
    rSig = x(ti+tau,1);
    tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
    tfr2(indices,icol)=rSig.*conj(dh(Lh+1+tau));
end

tfr1=fft(tfr1);
tfr2=fft(tfr2);

tfr1=tfr1(1:round(N/2),:);
tfr2=tfr2(1:round(N/2),:);

va=N/hlength;
IF=zeros(round(N/2),tcol);
% tfr=zeros(round(N/2),tcol);
E=mean(abs(x));

for i = 1 : round(N/2) % frequency
    for j = 1 : N      % time
        if abs(tfr1(i,j)) > 0.8*E%if you are interested in weak signals, you can delete this line.
            %if abs(1-real(va*1i*tfr2(i,j)/2/pi./tfr1(i,j)))<0.5
            if abs(-real(va*1i*tfr2(i,j)/2/pi./tfr1(i,j))) < 0.5
                IF(i,j) = 1;
            end
        end
    end
end
tfr = tfr1/(sum(h)/2);%the amplitude of tfr result has been pre-rectified.
Te = tfr.*IF;
% %The following code is an alternative way to estimate IF.
% %In theroy, they are same.
% omega = zeros(round(N/2),tcol);
% for b=1:N
% omega(:,b) = (ft-1)'+real(va*1i*tfr2(ft,b)/2/pi./tfr1(ft,b));
% end
% for i=1:round(N/2)%frequency
% for j=1:N%time
%     if abs(tfr1(i,j))>0.8*E%default frequency resolution is 1Hz.
%         if abs(omega(i,j)-i)<0.5%default frequency resolution is 1Hz.
%         IF(i,j)=1;
%         end
%     end
% end
% end

end

function [Ts] = SST_Y(x,hlength)
% Computes the SST (Ts)  of the signal x.
% INPUT
%    x      :  Signal needed to be column vector.
%    hlength:  The hlength of window function.
% OUTPUT
%    Ts     :  The SST
%   Written by YuGang in Shandong University at 2016.5.13.
[xrow,xcol] = size(x);

if (xcol~=1)
    error('X must be column vector');
end

if (nargin < 1)
    error('At least 1 parameter is required');
end

if (nargin < 2)
    hlength=round(xrow/5);
end

% Siglength=xrow;
hlength=hlength+1-rem(hlength,2);
ht = linspace(-0.5,0.5,hlength);ht=ht';

% Gaussian window
h = exp(-pi/0.32^2*ht.^2);
% derivative of window
dh = -2*pi/0.32^2*ht .* h; % g'

[hrow,~]=size(h); Lh=(hrow-1)/2;

N=xrow;
t=1:xrow;

[~,tcol] = size(t);


tfr1= zeros (N,tcol) ;
tfr2= zeros (N,tcol) ;

% tfr= zeros (round(N/2),tcol) ;
Ts= zeros (round(N/2),tcol) ;

for icol=1:tcol
    ti= t(icol); tau=-min([round(N/2)-1,Lh,ti-1]):min([round(N/2)-1,Lh,xrow-ti]);
    indices= rem(N+tau,N)+1;
    rSig = x(ti+tau,1);
    %rSig = hilbert(real(rSig));
    tfr1(indices,icol)=rSig.*conj(h(Lh+1+tau));
    tfr2(indices,icol)=rSig.*conj(dh(Lh+1+tau));
end

tfr1=fft(tfr1);
tfr2=fft(tfr2);

tfr1=tfr1(1:round(N/2),:);
tfr2=tfr2(1:round(N/2),:);

ft = 1:round(N/2);
bt = 1:N;

%%operator omega
nb = length(bt);
neta = length(ft);

va=N/hlength;
omega = zeros (round(N/2),tcol);

for b=1:nb
    omega(:,b) = (ft-1)'+real(va*1i*tfr2(ft,b)/2/pi./tfr1(ft,b));
end


omega=round(omega);

for b=1:nb%time
    % Reassignment step
    for eta=1:neta%frequency
        if abs(tfr1(eta,b))>0.0001%you can set much lower value than this.
            k = omega(eta,b);
            if k>=1 && k<=neta
                Ts(k,b) = Ts(k,b) + tfr1(eta,b);
            end
        end
    end
end
Ts=Ts/(xrow/2);
end
