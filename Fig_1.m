clear;clc;close all
Hz = 400; tend = 2;
t = 1/Hz : 1/Hz : tend;

Sig1 = exp(2*pi*1i*(-20*t.^2 + 90*t));               % IF1 = -40*t + 90;
Sig2 = exp(2*pi*1i*(20*t.^2 + 10*t));                % IF2 =  40*t + 10;          

Sig = Sig1 + Sig2;

% -------------------- STFT ----------------------------
window = 128; Nfrebin = 1024;
[Spec0,f] = STFT(real(Sig(:)),Hz,Nfrebin,window);
subplot(221);imagesc(t,f,abs(Spec0)/max(max(abs(Spec0)))); axis([0 2 0 100]);


xlabel('Time (s)'); ylabel('Frequency (Hz)');
set(gca,'YDir','normal'); cbar = colorbar; ylabel(cbar, 'Amplitude [•]'); 
caxis([0 1]);
set(cbar, 'Position', [0.5, 0.2, 0.02, 0.33]); 

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
xlabel('Time (s)'); ylabel('Number of components');
ylim([0 3])

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')