clc;close all; % clear;
% profile on
warning('off')
Hz = 400;
t = 0:1/Hz:2; t = t';

Sig1 = exp(2*pi*1i*(-20*t.^2 + 90*t));                   IF1 = -40*t + 90;
Sig2 = 0.6*exp(2*pi*1i*(20*t.^2 + 10*t));                IF2 =  40*t + 10;

Sig = Sig1 + Sig2;

Num = 2;    % the number of the components
deltaf = 20; % 0.03*Hz;
beta = 1e-3;

k = 0;
for Nc = 43 : 10 : 200
    j = 0; k = k + 1;
    for Nh = 46  : 10 : 200
        j = j + 1;

        % ridge extraction and fitting

        Sig0 = real(Sig)';
        if (isreal(Sig0))
            Sig0 = hilbert(Sig0);
        end

        fidexmult = zeros(Num,length(Sig0)); clearvars Spec Af c bw IFfit extr_Sig findex fidexmult
        num = Num;
        tic
        for i = 1 : Num
        
            [Spec,Atau,Af] = mALCT(Sig0,Hz,Nh,Nc,num);
            Spec = Spec';
            c = findridges(Spec,deltaf,0.5,num,7);

            c (c > round(length(t)/2)) = round(length(t)/2)-1;             % Added for rebuttal
            c (c < 1) = 1;                                                 % Added for rebuttal

            bw = Hz/60;           % 0.01*Hz;% the bandwidth of the TF filter for ridge extraction
            [IFfit,extr_Sig] = Dechirp_filter(Sig0,Hz,bw,Af(c),beta);
            fidexmult(i,:) = c;
            Sig0 = Sig0 - extr_Sig;
            num = num - 1;
        end
        toc

        % ------- ridge path regrouping

        Df = length(Af)/15; % Df = 0.03*Hz;
        [findex,interset] = RPRG(fidexmult,Df);

        % -------  Smoothing of ridges

        findex (findex < 1) = 1;                                           
        findex (findex > round(length(t)/2)) = round(length(t)/2)-1;       

        beta1 = 1e-5; % this parameter can be smaller which will be helpful for the convergence, but it may cannot properly track fast varying IFs
        iniIF = curvesmooth(Af(findex),beta1); % the initial guess for the IFs by ridge detection; the initial IFs should be smooth and thus we smooth the detected ridges

        e1 = norm(IF1-iniIF(1,:)')/norm(IF1);
        e2 = norm(IF2-iniIF(2,:)')/norm(IF2);

        e11 = norm(IF1-iniIF(2,:)')/norm(IF1);
        e21 = norm(IF2-iniIF(1,:)')/norm(IF2);

        e_k_1 = mean([e1 e2]);
        e_k_2 = mean([e11 e21]);

        e_k = min([e_k_1 e_k_2]);

        R80(k,j) = e_k;
        T(k,j) = toc;
    end
    k
end

disp(R80)
disp(T)

%%   Process data
A = zeros(161,161);
B = zeros(161,161);
%%
in = 4; m = 7; mm = m; hop = 48; detla_r = 15;
%%
rawi = 33; cols = 2;  rawf = rawi + detla_r;

for i = 1 : cols
    j = in; 
    for k = rawi: rawf
        A(j,m) = R36(k,i);
        j = j + 10;
    end
    m = m + 10;
end

rawi = rawi + hop; rawf = rawi + detla_r;

for i = 1 : cols
    j = in; 
    for k = rawi: rawf
        B(j,mm) = R36(k,i);
        j = j + 10;
    end
    mm = mm + 10;
end


%%  Plot Results

Nc = 40 : 200;
Nh = 40 : 200;
figure; set(gcf,'Position',[20 100 720 450]);
subplot(221);imagesc(Nh,Nc,A);
axis xy
strx = ('$N_h$');xlabel(strx,'Interpreter','LaTeX')
strx = ('$N_c$');ylabel(strx,'Interpreter','LaTeX')
cbar = colorbar; 
ylabel(cbar, 'IF error');
xlim([40 200]);ylim([40 200]);

subplot(222);imagesc(Nh,Nc,B);
axis xy
strx = ('$N_h$');xlabel(strx,'Interpreter','LaTeX')
strx = ('$N_c$');ylabel(strx,'Interpreter','LaTeX')
cbar = colorbar; 
ylabel(cbar, 'Time [s]'); 
xlim([40 200]);ylim([40 200]);

set(findall(gcf,'-property','FontSize'),'FontSize',10, 'FontName', 'Times New Roman')
%%
for nc1 = 1 : length(R80)
    A(10*(nc1-1)+9,2) = R80(nc1,1);
end
