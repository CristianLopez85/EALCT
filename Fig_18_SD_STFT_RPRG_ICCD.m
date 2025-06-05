clear;clc;close all

SampFreq = 400; tend = 2;
t = 1/SampFreq : 1/SampFreq : tend; Sf = 100;% time = time(1:400);
% 
amp1 = exp(0.1*t);  amp2 = exp(0.2*t);  amp3 = exp(0.3*t);
% 
Sig1 = amp1.*exp(2*1i*pi*(25*t.^2));                        IF1 = 50*t;
Sig2 = amp2.*exp(2*1i*pi*(25*t));                           IF2 = 25*ones(length(t),1);
Sig3 = amp3.*exp(2*1i*pi*(-22*t.^2 + 95*t + sin(5*t)));     IF3 = -44*t + 95 + 5*cos(5*t);

% Sig = Sig1 + Sig2;
%
Ratio = 0;
window = 128;
Nfrebin = 1024;
%
bw1 = SampFreq/60;%
orderIF1 = 50; % use the Fourier model to smooth the detected ridge curves£»orderIF1 could be larger than the following orderIF
num = 3; % the number of the components
delta = 0.03*SampFreq; %20; % ;                                                            % attempt by CL 0.03*SampFreq
alpha1 = 5;
%
alpha2 = 0.5;
bw = SampFreq/40; % 2.5; %  % bandwidth of the ICCD                          
orderIF = 5;
%
% beta1 = 1e-1; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs

j = 0; pin = 0; % 200 to compare with EALT
for SNR = -4:1:15   % do separately
    j = j + 1;
    for k = 1: 10
        pin = pin + 1;
        Sig = Sig1 + Sig2 + Sig3;
        seed = pin; % Use the loop variable as the seed
        rng(seed); % Set the seed for random number generation
        Sig = awgn(Sig,SNR);

        %%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
        [Spec,f] = STFT(Sig(:),SampFreq,Nfrebin,window);

        % ridge extraction and fitting
        [fidexmult, tfdv] = extridge_mult(Sig, SampFreq, num, delta, orderIF1,bw1,Nfrebin,window,alpha1);

        % ridge path regrouping
        thrf = 0.03*SampFreq; % length(f)/30; %                                             %  attempt by by CL  0.03*SampFreq 
        
        [findex,interset] = RPRG(fidexmult,thrf);

%         findex1 = round(curvesmooth(findex,beta1)); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges
% 
%         figure;subplot(121);plot(fidexmult(2,:));hold on; plot(findex(1,:),'r--');plot(findex1(1,:));legend('ridge','IF','IF-Smooth')
%         subplot(122);plot(fidexmult(2,:));hold on; plot(findex(2,:),'r--');plot(findex1(2,:));legend('ridge','IF','IF-Smooth')

        % signal decomposition
        orderamp = round(bw*length(Sig)/SampFreq);% Fourier order for characterizing signal amplitudes
        [extr_Sig,ampmatrix,IFfit] = ICCD(Sig,SampFreq,f(findex),orderIF,orderamp,alpha2); % ICCD performs the joint-optimization scheme using all the obtained IFs
        
        % by CL
        chops = 100;
        Sig1_s = Sig1(chops: end-chops); Sig2_s = Sig2(chops: end-chops); Sig3_s = Sig3(chops: end-chops);
        mode_s = extr_Sig(1,chops: end-chops); mode2_s = extr_Sig(2,chops: end-chops); mode3_s = extr_Sig(2,chops: end-chops);

        e1 = norm(real(Sig1_s-mode_s))/norm(real(mode_s));
        e2 = norm(real(Sig2_s-mode2_s))/norm(real(mode2_s));
        e3 = norm(real(Sig3_s-mode3_s))/norm(real(mode3_s));

        e1_r = norm(real(Sig2_s-mode_s))/norm(real(mode_s));               % Because, may Sig1_s can be displayed as Sig2_s
        e2_r = norm(real(Sig3_s-mode2_s))/norm(real(mode2_s));
        e3_r = norm(real(Sig1_s-mode3_s))/norm(real(mode3_s));

        e1_r1 = norm(real(Sig3_s-mode_s))/norm(real(mode_s));               % Because, may Sig1_s can be displayed as Sig2_s
        e2_r1 = norm(real(Sig1_s-mode2_s))/norm(real(mode2_s));
        e3_r1 = norm(real(Sig2_s-mode3_s))/norm(real(mode3_s));
        
        e_M1 = [e1 e1_r e1_r1];
        e_M2 = [e2 e2_r e2_r1];
        e_M3 = [e3 e3_r e3_r1];

        e1_s = min(e_M1); e2_s = min(e_M2); e3_s = min(e_M3);

        e_k(k) = mean([e1_s e2_s e3_s]);
    end
    e(j) = mean(e_k);
    j
end
e=e';

%%
figure; 
set(gcf,'Position',[20 100 800 450]);	
subplot(221);plot(t,real(Sig1),'k');hold on;plot(t,real(extr_Sig(1,:)),'--c');ylim([-1.5 2])
xlabel('Time (s)');ylabel('Amplitude');legend('Original','Reconstructed','Orientation','horizontal');
subplot(223);plot(t,real(Sig1-extr_Sig(1,:)),'r');ylim([-1.5 2]);xlabel('Time (s)');ylabel('Amplitude');legend('Error')
subplot(222);plot(t,real(Sig2),'k');hold on;plot(t,real(extr_Sig(2,:)),'--c');ylim([-1.5 2]);
xlabel('Time (s)');ylabel('Amplitude');legend('Original','Reconstructed','Orientation','horizontal');
subplot(224);plot(t,real(Sig2-extr_Sig(2,:)),'r');xlabel('Time (s)');ylabel('Amplitude');
xlabel('Time (s)');ylabel('Amplitude');legend('Error');ylim([-1.5 2]);

set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')