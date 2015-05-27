%close all
clear all;


%% system to be controlled

mode = 'azuma08fwd';           % Azuma 08 Automatica
%mode = 'azuma08feedback';      % uses Azuma 08 Automatica in feedback

%mode = 'azumaLPfwd';           % uses Azuma 08 TAC finite-horizon LP
%mode = 'azumaLPfeedback';      % uses Azuma 08 TAC in feedback

%mode = 'minami07';             % uses Minami 07 CDC (specific to feedback)


Tvec = 30;

maxBound = zeros(1,length(Tvec));
for j = 1:length(Tvec)
    
    if(length(Tvec)>3)
        sysParams.plot=0;
    else
        sysParams.plot=1;
    end
    if( strcmp(mode,'azumaLPfwd') || strcmp(mode,'azumaLPfeedback') )
        %sysParams.T = 10;        % horizon length for LP
        sysParams.T = Tvec(j);
        %sysParams.name = 'automatica08ex';
        sysParams.gamma = 100;    % bound on norm of Q_wv
    else
        sysParams.T = [];
        sysParams.gamma = [];
    end
    % QUANTIZER INTERVAL
    d = 2;
    
    % system - NOTES:
    % Bp (discrete B for input) is B2 in Minami 07
    % Cp (discrete y = Cx) is C1 in Minami 07
    % State feedback gain (F) is C2 in Minami 07
    
    %{
    % system from Azuma 08 Automatica:
    Apc = [-3 3;0 -2];
    Bpc = [0;2];
    Cpc = [1 1];
    Dpc = 0;
    x0 = [1;2];
    %}
    
    %
    % simple m, k, b sys
    mm = 1;kk = 1;bb = .1;
    Apc = [0 1;(-kk/mm) (-bb/mm)];
    Bpc = [0;1/mm];
    %%Cpc = [1 0];
    Cpc = eye(2);
    Dpc = 0;
    x0 = [1;2];
    %}
    
    
    % convert to discrete time
    dt = 0.1;
    Tf = 15; N = ceil(Tf/dt)+1;
    %N = 100;
    
    sysC = ss(Apc,Bpc,Cpc,Dpc);
    sysParams.x0 = x0;
    sysParams.dt = dt;
    sysParams.N = N;
    sysParams.d = d;
    
    if( strcmp(mode,'azuma08feedback') || strcmp(mode,'azumaLPfeedback') )
        % state feedback gains: (design these based on robust control and w?)
        F = [0.007 0.859];
        sysParams.F = F;
        % observer-based controller (for 08 Automatica) gains:
        L = [-0.032 0.202]';
        sysParams.L = L;
    end
    
    
    %{

%%%% TESTING QUANTIZER

[Qstar,OPT] = designAzuma08TAC(sysC,sysParams);

T = sysParams.T;
Aq = Qstar.A;
Bq2 = Qstar.B2;
Cq = Qstar.C;

% Bg = Qstar.Bg;
% Cg = Qstar.Cg;
% Ag = Qstar.Ag - Bg*Cg;

% Test impulse response of system realization
for ktest = 0:(T-1)
    std = (Cq*(Aq + Bq2*Cq)^ktest*Bq2);
%     Himp2ss = (Cg*(Ag + Bg*Cg)^ktest*Bg);
    H2k = OPT.H2(ktest+1);
    fprintf('std: %f   H2k: %f\n',std,H2k)
end
    
    %}
    
    
    
    
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    o = azuma_test(mode,sysC,sysParams);

    %{
    fprintf('\nAq:\n')
    disp(o.Qstar.A)
    fprintf('B2q\n')
    disp(o.Qstar.B2)
    fprintf('Cq:\n')
    disp(o.Qstar.C)
    %}
    
    %
% (for convenience)
Pd = c2d(sysC,dt);[Ap,Bp,Cp,~] = ssdata(Pd);[~,m] = size(Bp);[p,~] = size(Cp);
T = sysParams.T;
[n,~] = size(Ap);
nQ = o.Qstar.N;Cq = o.Qstar.C;Bq2 = o.Qstar.B2;Aq = o.Qstar.A;
Cbar = [Cp zeros(p,nQ)];
Abar = [Ap,Bp*Cq;zeros(nQ,n),Aq+Bq2*Cq];
B2bar = [Bp;Bq2];
Phi = zeros(p,m,T);
for k = 1:T
    Phi(:,:,T-k+1) = Cp*Ap^(k-1)*Bp;
end


%
% check bound with bar, opt vars
sumCABbar = zeros(p,T-1);
sumEquiv = zeros(p,T-1);
for k=1:T-1;
    sumCABbar(:,k) = (Cbar*Abar^k*B2bar);
    CAB = Cp*Ap^k*Bp;
    PhiH2 = sum(squeeze(Phi(:,:,(T-k+1):T))*squeeze(o.OPTresults.H2(:,:,1:k)));
    sumEquiv(:,k) = CAB + PhiH2;
end
%sum(abs(sumCAB))
%sum(abs(sumEquiv))

maxBound(j) = norm( abs(Cp*Bp) + sum(abs(sumCABbar),2) ,'inf')*(d/2);
%CAB_bar = Cbar*Abar^k*B2bar
%equiv = Cp*Ap^k*Bp + sum(Phi(:,:,(T-k+1):T).*o.OPTresults.H2(:,:,1:k))

%}
    
    
end




figure
plot(Tvec,maxBound)
xlabel('horizon T')
ylabel('E(Q*)')
title('Bound on input-output quantization error')
