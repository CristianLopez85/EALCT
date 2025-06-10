clear;clc;close all
tic
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


%%
j = 0;pin = 0;
for SNR = -4:15
    j = j + 1;
    for k = 1:10
        pin = pin + 1;

        clearvars -except Hz j k Sig1 Sig2 Sig3 SNR t e_k e j pin
        seed = pin; % Use the loop variable as the seed
        rng(seed); % Set the seed for random number generation

        Sig = Sig1 + Sig2 + Sig3;
%         Sig = awgn(Sig,SNR);

        beta = 1e-3;
        Sig0 = real(Sig)';
        if (isreal(Sig0))
            Sig0 = hilbert(Sig0);
        end
        Num = 3;                % the number of the components
        fidexmult = zeros(Num,length(Sig0)); clearvars Spec Af c bw IFfit extr_Sig findex fidexmult
        Nh = 151; Nc = 71;
        num = Num;
        for i = 1 : Num

            [Spec,Atau,Af] = zALCT(Sig0,Hz,Nh,Nc,num);
            Spec = Spec';
            deltaf = 0.03*Hz; % 20;            % before: ; 20 in Chen's code
            c = findridges_m(Spec,deltaf,0.5,num,7);

            c (c > round(length(t)/2)) = round(length(t)/2)-1;             % Added for rebuttal
            c (c < 1) = 1;                                                 % Added for rebuttal

            % ridge extraction and fitting
            bw = 0.01*Hz; % Hz/60; % ; ; Hz/60 in Chen's code
            [IFfit,extr_Sig] = Dechirp_filter(Sig0,Hz,bw,Af(c),beta);
            fidexmult(i,:) = c;
            Sig0 = Sig0 - extr_Sig;
            num = num - 1;
        end

        % ridge path regrouping
        Df = 0.03*Hz; % length(Af)/15; % ;  length(f): in Chen's code
        [findex,interset] = RPRG(fidexmult,Df);

        findex (findex < 1) = 1;                                           % Added for rebuttal
        findex (findex > round(length(t)/2)) = round(length(t)/2)-1;       % Added for rebuttal


        % Smoothing of ridges
        beta1 = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
        iniIF = curvesmooth(Af(findex),beta1); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges

        % run the CT with g for TFC domain
        alpha1 = 50;
        tt = -8:1/Hz:8;
        h1 = exp(-pi*alpha1*tt.^2'); % window g_0
        eta = 1;
        [tfc1, tfrtic, tcrtic] = sqSTCTm(Sig, eta/length(Sig), 1, h1);

        % TFC ridges
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

        line1 (line1 < 1) = 1;                                                 % Added for rebuttal
        line2 (line2 < 1) = 1;                                                 % Added for rebuttal
        line3 (line3 < 1) = 1;                                                 % Added for rebuttal


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

        % Rebuttal

        chops = 100;
        Sig1_s = Sig1(chops: end-chops); Sig2_s = Sig2(chops: end-chops); Sig3_s = Sig3(chops: end-chops);
        mode = rrcon(1,:)'; mode2 = rrcon(2,:)'; mode3 = rrcon(3,:)';
        mode_s = mode(chops: end-chops); mode2_s = mode2(chops: end-chops); mode3_s = mode3(chops: end-chops);


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
        % k
    end
        e(j) = mean(e_k);
        j
end
% figure;
% subplot(311);plot(real(Sig1_s));hold on;plot(real(mode_s),'r--')
% subplot(312);plot(real(Sig2_s));hold on;plot(real(mode2_s),'r--')
% subplot(313);plot(real(Sig3_s));hold on;plot(real(mode3_s),'r--')
% 
% figure;
% subplot(311);plot(real(Sig2_s));hold on;plot(real(mode_s),'r--')
% subplot(312);plot(real(Sig3_s));hold on;plot(real(mode2_s),'r--')
% subplot(313);plot(real(Sig1_s));hold on;plot(real(mode3_s),'r--')
% 
% figure;
% subplot(311);plot(real(Sig3_s));hold on;plot(real(mode_s),'r--')
% subplot(312);plot(real(Sig1_s));hold on;plot(real(mode2_s),'r--')
% subplot(313);plot(real(Sig2_s));hold on;plot(real(mode3_s),'r--')

%%
% figure; set(gcf,'Position',[20 100 800 450]);
% subplot(231)
% plot(t,real(Sig1), 'color', 'k'); xlabel('time (s)');hold on
% plot(t,real(rrcon(3,:)),'--c');ylim([-2.1 2.1]);
% xlabel('Time (s)'); ylabel('Amplitude');%legend('Real','Reconstructed','Orientation','horizontal');
% 
% subplot(234)
% plot(t, real(Sig1')-real(rrcon(3,:)), 'r'); ylim([-2.1 2.1]); xlabel('Time (s)');ylabel('Amplitude');legend('Error')
% 
% subplot(232)
% plot(t,real(Sig2), 'color', 'k'); xlabel('Time (s)');hold on;
% plot(t,real(rrcon(2,:)),'--c');ylim([-2.1 2.1]);
% xlabel('Time (s)'); ylabel('Amplitude');legend('Real','Reconstructed','Orientation','horizontal');
% 
% subplot(235)
% plot(t, real(Sig2')-real(rrcon(2,:)), 'r'); ylim([-2.1 2.1]); xlabel('Time (s)');ylabel('Amplitude');legend('Error')
% 
% subplot(233)
% plot(t,real(Sig3), 'color', 'k'); xlabel('Time (s)');hold on;
% plot(t,real(rrcon(1,:)),'--c');ylim([-2.1 2.1]);
% xlabel('Time (s)'); ylabel('Amplitude');%legend('Real','Reconstructed','Orientation','horizontal');
% 
% subplot(236)
% plot(t, real(Sig3')-real(rrcon(1,:)), 'r'); ylim([-2.1 2.1]); xlabel('Time (s)');ylabel('Amplitude');legend('Error')
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
% figure;plot(e)
toc
e = e';
