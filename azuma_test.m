function output = azuma_test(mode,sysC,sysParams)
% runs test
% BR, 2/13/2013

x0 = sysParams.x0;
dt = sysParams.dt;
N = sysParams.N;
d = sysParams.d;    % quantizer interval:
if(~isfield(sysParams,'plot'))
    sysParams.plot=1;
end

% convert to discrete time
Tf = (N-1)*dt;
t = 0:dt:Tf;
Pd = c2d(sysC,dt);[Ap,Bp,Cp,~] = ssdata(Pd);
[n,~] = size(Ap);[p,~] = size(Cp);[~,m] = size(Bp);

if( strcmp(mode,'azuma08feedback') || strcmp(mode,'azumaLPfeedback') )
    F = sysParams.F;    % state feedback gain
    L = sysParams.L;
    % estimator/observer gains:
    Ao = (Ap - L*Cp);
    xhat = x0;
end

if( (strcmp(mode,'azuma08fwd')) || (strcmp(mode,'azuma08feedback')) )
    % CHECK ASSUMPTIONS
    if(p~=m);disp('ERROR: input-output dimensions not consistent');end
    if(eig(Cp*Bp)==0);disp('ERROR: (CB) is singular');end
    
    % design quantizer (Based on Azuma, 2008 Automatica)
    Aq = Ap;    Bq1 = -Bp;    Bq2 = -Bq1;
    Cq = -inv(Cp*Bp)*Cp*Ap;    Dq = 1;
    nQ = n;
    Qstar.A = Aq;    Qstar.B1 = Bq1;    Qstar.B2 = Bq2;
    Qstar.C = Cq;    Qstar.N = nQ;
    

elseif( strcmp(mode,'minami07'))
    %find tau:
    tau = 0;while( Cp*(Ap+Bp*F)^tau*Bp == 0);tau = tau+1;end
    Aq = Ap + Bp*F;    Bq1 = -Bp;    Bq2 = -Bq1;
    Cq = -pinv(Cp*Aq^tau*Bp)*Cp*Aq^(tau+1);    Dq = 1;
    nQ = n;
    Qstar.A = Aq;    Qstar.B1 = Bq1;    Qstar.B2 = Bq2;
    Qstar.C = Cq;    Qstar.N = nQ;
    
    % CHECK ASSUMPTIONS
    if( abs(eig(Ap + Bp*F)) > 1);disp('Error - unstable closed loop');end
    if(rank(F)<q);disp('ERROR: F (C2) Not full row rank');end
    if( rank(Cp*(Ap+Bp*F)^tau*Bp) < p )
        fprintf('\nERROR: Cp*(Ap+Bp*F)^tau*Bp not full row rank, tau = %d\n\n',tau)
    end
    
elseif( strcmp(mode,'azumaLPfwd') || strcmp(mode,'azumaLPfeedback') )
    [Qstar,OPTresults] = designAzuma08TAC(sysC,sysParams);
    nQ = Qstar.N;
    Aq = Qstar.A;
    Bq1 = Qstar.B1;
    Bq2 = Qstar.B2;
    Cq = Qstar.C;
    Dq = 1;
    output.OPTresults = OPTresults;
    
end
output.Qstar = Qstar;

%% simulation

% preallocate based on sim length
xSave = zeros(n,N);     % state: dyanamically quantized control (v)
xcSave = zeros(n,N);    % state: continuous control input (u)
xUSave = zeros(n,N);    % state: naive control quantization (vU)
xiSave = zeros(nQ,N);    % dynamic quantizer state
uSave = zeros(m,N);     % continuous control
vSave = zeros(m,N);     % dynamically quantized control
vUSave = zeros(m,N);    % naive quantized control
xHatSave = zeros(n,N);  % (azuma 08 with feedback): observer state

% initial conditions
x = x0;xU = x0;xc = x0;
xSave(:,1) = x0;xUSave(:,1) = x0;xcSave(:,1) = x0;
xi = zeros(nQ,1);

for k = 1:(N-1)
    
    switch mode
        
        case 'azuma08fwd'       % FEEDFORWARD CONTROLLER
            u = sin(0.1*pi*k) + 2*cos(0.07*pi*k);
            argq = Cq*xi + Dq*u;
            
        case 'azumaLPfwd'
            u = sin(0.1*pi*k) + 2*cos(0.07*pi*k);
            argq = Cq*xi + Dq*u;
            
        case 'azuma08feedback'  % FEEDBACK CONTROLLER, WITH OBSERVER
            
            if(k==1);  u=0;  end
            % run estimator/observer
            y = Cp*x;
            xhat = Ao*xhat + L*y + (Bp-L*Cp*Bp)*u;
            
            % compute continuous control
            u = F*xhat;
            argq = Cq*xi + Dq*u;
            
            % save vars
            xHatSave(:,k) = xhat;
            
        case 'azumaLPfeedback'
            if(k==1);  u=0;  end
            % run estimator/observer
            y = Cp*x;
            xhat = Ao*xhat + L*y + (Bp-L*Cp*Bp)*u;
            
            % compute continuous control
            u = F*xhat;
            argq = Cq*xi + Dq*u;
            
            % save vars
            xHatSave(:,k) = xhat;
            
        case 'minami07'
            
            u = F*x;
            argq = Cq*xi + u;
            
    end
    
    % quantize control
    v = staticNearestNeighbor(argq,d);
    vU = staticNearestNeighbor(u,d);
    
    % apply quantized control, update true state (with both inputs)
    x = Ap*x + Bp*v;
    xc = Ap*xc +Bp*u;
    xU = Ap*xU + Bp*vU;
    
    % update dynamic quantizer state
    xi = Aq*xi + Bq1*u + Bq2*v;
    
    uSave(k) = u;    vSave(k) = v;    vUSave(k) = vU;
    xSave(:,k+1) = x;    xUSave(:,k+1) = xU;    xcSave(:,k+1) = xc;
    xiSave(:,k+1) = xi;
    
end

%%
if( (strcmp(mode,'azuma08fwd')) || (strcmp(mode,'azuma08feedback')) )
    maxBound = max(abs(Cp*Bp))*d/2;
elseif( (strcmp(mode,'minami07')) )
    maxBound = max(abs(Cp*Aq^tau*Bp))*d/2;
elseif( (strcmp(mode,'azumaLPfwd')) || (strcmp(mode,'azumaLPfeedback')) )
    %try
        Cbar = [Cp zeros(p,nQ)];
        Abar = [Ap,Bp*Cq;zeros(nQ,n),Aq+Bq2*Cq];
        B2bar = [Bp;Bq2];
        
        sumCAB = 0;
        for k=1:(sysParams.T-1);
            sumCAB = sumCAB + abs(Cbar*Abar^k*B2bar);
        end
        maxBound = norm( abs(Cp*Bp) + sumCAB ,'inf')*(d/2);
        fprintf('\nmax bound E(Q*) = %f \n\n',maxBound)
%     catch
%         disp('LP bound bug')
%         maxBound = 0;
%     end
end

norm2v = norm( (Cp*xcSave - Cp*xSave), 2);
norm2u = norm( (Cp*xcSave - Cp*xUSave), 2);
normInfv = norm( (Cp*xcSave - Cp*xSave), 'inf');
normInfu = norm( (Cp*xcSave - Cp*xUSave), 'inf');

output.maxBound = maxBound;
output.norm2v = norm2v;
output.norm2u = norm2u;
output.normInfv = normInfv;
output.normInfu = normInfu;

if(sysParams.plot)
figure
subplot(3,1,1)
stairs(t,uSave,'b')
hold on
stairs(t,vUSave,'g')
stairs(t,vSave,'r')
if( strcmp(mode,'azumaLPfwd') || strcmp(mode,'azumaLPfeedback') )
    title(sprintf('Control, T = %d',sysParams.T))
else
    title('Control')
end

subplot(3,1,2)
h0 = stairs(t,(Cp*xcSave)','b');
hold on
h1 = stairs(t,(Cp*xUSave)','g');
h2 = stairs(t,(Cp*xSave)','r');
title('Output')
legend([h0(1) h1(1) h2(1)],'y_{u}','y_{q(u)}','y_{v}')

subplot(3,1,3)
h0 = plot([min(t) max(t)],[maxBound maxBound],'k--');
hold on
plot([min(t) max(t)],[-maxBound -maxBound],'k--')
h1 = stairs(t,(Cp*xcSave - Cp*xUSave)','g');
h2 = stairs(t,(Cp*xcSave - Cp*xSave)','r');
title(sprintf('  L2 error:     q(u): %0.2f,  v: %0.2f \nMax error:   q(u): %0.2f,  v: %0.2f',norm2u,norm2v,normInfu,normInfv))
legend([h1(1) h2(1) h0],'y(u)-y(q(u))','y(u)-y(v)','E(Q*)')
xlabel('time [s]')
end

%figure
%stairs(xiSave')


end





