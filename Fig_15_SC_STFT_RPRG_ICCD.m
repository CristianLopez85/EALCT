clc; clear; close all
SampFreq = 100;
t = 0:1/SampFreq:3; 

Sig1 = exp(2*pi*1i*(30*t + 3*sin(3*t)));                                   %  IF2 = 30 + 9*cos(3*t);
Sig2 = 0.9*exp(2*pi*1i*(30*t - 3*sin(3*t) + 3/(2*pi)));                    %  IF1 = 30 - 9*cos(3*t);

% Sig = Sig1 + Sig2;
%
Ratio = 0;
window = 128/2;
Nfrebin = 1024/2;
%
bw1 = SampFreq/60;%                                          % attempt by CL 0.01*SampFreq
orderIF1 = 50; % use the Fourier model to smooth the detected ridge curves£»orderIF1 could be larger than the following orderIF
num = 2; % the number of the components
delta = 0.03*SampFreq;                                                      % attempt by CL 0.03*SampFreq
alpha1 = 5;
%
alpha2 = 0.5;
bw = SampFreq/40; % 2.5; %  % bandwidth of the ICCD                         
orderIF = 5;
%
% beta1 = 1e-1; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs

j = 0; pin = 0; 
for SNR = -4:1:15   % do separately
    j = j + 1;
    for k = 1: 10
        pin = pin + 1;
        Sig = Sig1 + Sig2;
        seed = pin; % Use the loop variable as the seed
        rng(seed); % Set the seed for random number generation
        Sig = awgn(Sig,SNR);

        %%%%%%%%%%%STFT%%%%%%%%%%%%%%%%%%%%%
        [Spec,f] = STFT(Sig(:),SampFreq,Nfrebin,window);

        % ridge extraction and fitting
        [fidexmult, tfdv] = extridge_mult(Sig, SampFreq, num, delta, orderIF1,bw1,Nfrebin,window,alpha1);

        % ridge path regrouping
        thrf = 0.03*SampFreq;%0.1*SampFreq; %length(f)/30; %                                             %  attempt by by CL  0.03*SampFreq 
        
        [findex,interset] = RPRG(fidexmult,thrf);

%         findex1 = round(curvesmooth(findex,beta1)); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges
% 
%         figure;subplot(121);plot(fidexmult(1,:));hold on; plot(findex(1,:),'r--');plot(findex1(1,:));legend('ridge','IF','IF-Smooth')
%         subplot(122);plot(fidexmult(2,:));hold on; plot(findex(2,:),'r--');plot(findex1(2,:));legend('ridge','IF','IF-Smooth')

        % signal decomposition
        orderamp = round(bw*length(Sig)/SampFreq);% Fourier order for characterizing signal amplitudes
        [extr_Sig,ampmatrix,IFfit] = ICCD(Sig,SampFreq,f(findex),orderIF,orderamp,alpha2); % ICCD performs the joint-optimization scheme using all the obtained IFs

        % by CL
        chop = 25;
        Sig1_s = Sig1(chop: end-chop); Sig2_s = Sig2(chop: end-chop);
        mode_s = extr_Sig(1,chop: end-chop); mode2_s = extr_Sig(2,chop: end-chop);

        e1 = norm(real(Sig1_s-mode_s))/norm(real(mode_s));
        e2 = norm(real(Sig2_s-mode2_s))/norm(real(mode2_s));

        e1_r = norm(real(Sig2_s-mode_s))/norm(real(mode_s));               % Because, may Sig1_s can be displayed as Sig2_s
        e2_r = norm(real(Sig1_s-mode2_s))/norm(real(mode2_s));

        e_k_1 = mean([e1 e2]);
        e_k_2 = mean([e1_r e2_r]);

        e_k(k) = min([e_k_1 e_k_2]);
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