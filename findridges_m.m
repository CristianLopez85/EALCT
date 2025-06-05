function [index,limits] = findridges_m(Spec,delta,penalty,k,Np)
%Ridge detection algorithm,i.e., Algorithm 1 in paper:Separation of Overlapped Non-Stationary Signals by Ridge Path Regrouping and Intrinsic Chirp Component Decomposition
% ,IEEE Sensors journal,2017.
% The algorithm is originally introduced in the paper: Algorithms for blind components separation and extraction from the time-frequency distribution of their mixture.
% EURASIP Journal on Advances in Signal Processing, 2004.


%Spec£ºTime-Frequency distribution of the signal
%delta£ºmaximum allowable frequency variation between two consecutive points
%index£ºThe obtained frequency indexs at each time instant

Spec = abs(Spec);
[M,N] = size(Spec);
% index = zeros(1,N);
index = k*ones(1,N);                          % --------------- Added by CL
[fmax,tmax] = find(Spec == max(Spec(:)));
fmax = fmax(1);
tmax = tmax(1);
index(tmax) = fmax;

f0 = fmax;     

limits = zeros(1,2);
E = sum(Spec); mu = min(E)/max(E);       % ---- Added by CL, from Miao
D_E = penalty*mu*max(max(Spec));        % ---- Added by CL, from Miao
for j = (min(tmax+1,N)):N
    low = max(1,f0-delta);
    up = min(M,f0+delta);
    [Am,f0] = max(Spec(low:up,j));
    f0 = f0 + low - 1;
    index(j) = f0;
    %                                         
    if Am < D_E                              % ---- Added by CL, from Miao
        b = index(1:j-1);                    % --------------- Added by CL
        m = (b(j-1)-b(j-1-Np))/Np;   %
        x = linspace(j,N,N-(j)+1);
        y = round(m*(x-j-1)+b(j-1));             % linear fitting
        index = [b(1:j-1) y];
%         j
        break
    end
    % ---------------
end
limits(1,2) = j;                              % --------------- Added by CL

f1 = fmax;

for j = (max(1,tmax-1)):-1:1
    low = max(1,f1-delta);
    up = min(M,f1+delta);
    [Am,f1] = max(Spec(low:up,j));
    f1 = f1 + low - 1;
    index(j) = f1;
    %                                         % --------------- Added by CL
    if Am < D_E
        b = index(j+1:N);
        m = (b(1+Np)-b(1))/Np;
        x = linspace(1,j,j);
        y = round(m*(x-j+1)+b(1));             % linear fitting
        index = [y b];
        break
    end
    % ---------------
end
% index(index<k)=k;
limits(1,1) = j;                              % --------------- Added by CL
end
