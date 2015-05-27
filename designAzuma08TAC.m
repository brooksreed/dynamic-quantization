function [Qstar,OPTresults] = designAzuma08TAC(sysC,sysParams)
% designs dynamic quantizer using Azuma 08 TAC (LP) method
% function [Qstar,OPTresults] = designAzuma08TAC(sysC,sysParams)
% sysC: continuous time state space syste
% sysParams: struct with dt, T (horizon)
% Qstar: struct with (N,A,B1,B2,C) optimal quantizer 
% OPTresults: LP solution (impulse response matrices H2k)

% BR, 3/4/2013
% changelog: author,date,change
%{
- 3/6/2013: fixed bug in SVD for Ho-Kalman (Wc')
-

%}

T = sysParams.T;dt = sysParams.dt;
Pd = c2d(sysC,dt);[~,Bp,~,~] = ssdata(Pd);
%[p,~] = size(Cp);  % p OUTPUTS
[~,m] = size(Bp);   % m INPUTS

% solve with AMPL
OPTresults = solveAMPL_AzumaLP(sysC,sysParams);

% formulate nominal quantizer based on (43)-(46) 
% first construct block Hankel H
TT = 2*(T-1);allH = zeros(m,m,TT);
allH(:,:,1:T) = OPTresults.H2;
allH(:,:,T+1:TT) = zeros(m,m,TT-T);

% allH = OPTresults.H2;

[H,Tp] = makeBlockHankel(allH);
% svd
[Wo,S,Wc] = svd(H);
Wc = Wc';   % U*S*V'
OPTresults.H = H;OPTresults.Wo = Wo;OPTresults.S = S;OPTresults.Wc = Wc;

% formulate quantizer
No = m*Tp;
B2o = (S)^(1/2)*Wc*[eye(m)' zeros(m,m*(Tp-1))]';
Co = [eye(m) zeros(m,No-m)]*Wo*(S)^(1/2);
Ao = pinv([eye(m*(Tp-1)) zeros(m*(Tp-1),m)]*Wo*(S)^(1/2))*...
    [zeros(m*(Tp-1),m) eye(m*(Tp-1))]*Wo*(S)^(1/2) - B2o*Co;

% check stability of Ao + B2o*Co
if(max(abs(eig(Ao + B2o*Co)))<1)
    % QSTAR structure
    Qstar.N = No;
    Qstar.A = Ao;
    Qstar.B1 = -B2o;
    Qstar.B2 = B2o;
    Qstar.C = Co;
    disp('nominal quantizer stable - reduced form')
    %disp(eig(Ao+B2o*Co))
    
else % nominal quantizer unstable - need augmented quantizer
    disp('nominal quantizer unstable - augmented form')
    %disp(eig(Ao+B2o*Co))
    
    Nobar = No*T;
    B2obar = [repmat(zeros(size(B2o')),1,(T-1)) B2o']';
    Cobar = repmat(Co,1,T);
    Aobar = zeros(Nobar);
    [Borows,~] = size(B2o);
    
    Aobar((Nobar-Borows+1):Nobar,:) = repmat(-(B2o*Co),1,T);
    Aobar(1:(Nobar-Borows),:) = [zeros(Nobar-Borows,Borows) kron(eye(T-1),Ao-B2o*Co)];
    
    %eig(Aobar)
    %eig(Aobar + B2obar*Cobar)
    
    Qstar.N = Nobar;
    Qstar.A = Aobar;
    Qstar.B1 = -B2obar;
    Qstar.B2 = B2obar;
    Qstar.C = Cobar;
    
    
    %%%
%     for ktest = 0:(T-1)
%         fprintf('\nbar:\n')
%         disp(Cobar*(Aobar + B2obar*Cobar)^ktest*B2obar)
%         fprintf('std:\n')
%         disp(Co*(Ao + B2o*Co)^ktest*B2o)
%         fprintf('H2k:\n')
%         disp(OPTresults.H2(ktest+1))
%     end

end

end



function results = solveAMPL_AzumaLP(sysC,sysParams)

% BR, 3/4/2013
% changelog: author,date,change
%{
-
%}

% solve options
%generalOptions (all defaults = 0)
%   struct with
%   '.noAMPLPresolve','.noCPLEXPresolve','.AMPLtimes'
%   '.parallelmode','.threads','.timeoutCPLEX'
generalOptions.parallelmode=-1; % -1 opportunistic, 1 deterministic
generalOptions.threads=0;
generalOptions.timeoutCPLEX=3*60;
solveOptions=[];

% default solve options for LP
if(isempty(solveOptions))
    solveOptions.alg='dual';
    solveOptions.PD='primal';
end

%{
% cd to matlab MET folder
resultsFolder=...
    'C:\Brooks\Dropbox\Research Dropbox\MATLAB Code\Quantization\LP Results';
current = cd(resultsFolder);

% start diary
diary off;
diary('cmdLog.txt')

% prep save files
if(isfield(sysParams,'name'))
    fPrefix=[sysParams.name, '_']; else fPrefix = 'test_'; end

% save filename for problem instance (empty string if no save)
saveFilename=sprintf('%sT%i_%s',fPrefix,sysParams.T,...
    dateString('DHMS'));
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cd to location of AMPL models, temp data and saved data/outputs
current = cd...
    ('C:\Brooks\Dropbox\Research Dropbox\Azuma LP Raw Files');

% generate .dat
[p,m,T] = generateData_AzumaLP(sysC,sysParams);

%make ampl .run file
fid = fopen('AzumaLP.run', 'w');

writeAMPLRunOptions(fid,'LP',generalOptions,solveOptions);
fprintf(fid, 'model AzumaLP.mod;\n');
fprintf(fid, 'data AzumaLP.dat;\n');

% SOLVE with CPLEX
fprintf(fid, 'solve;\n\n');

% AMPL run file: printfs for data to be read-in by matlab: H2,H2bar,Epsbar

% H2
fprintf(fid, 'for {i in 1..m} {\n');
fprintf(fid, '   for {j in 1..m} {\n');
fprintf(fid, '       for {k in 1..T} {\n');
fprintf(fid,...
    '            printf " %%i,%%i,%%i,%%f\\n",i,j,k,H2[i,j,k] > H2Out.ans;\n');
fprintf(fid, '       }\n');
fprintf(fid, '   }\n');
fprintf(fid, '}\n');

% H2bar
fprintf(fid, 'for {i in 1..m} {\n');
fprintf(fid, '   for {j in 1..m} {\n');
fprintf(fid, '       for {k in 1..T} {\n');
fprintf(fid,...
    '            printf " %%i,%%i,%%i,%%f\\n",i,j,k,H2bar[i,j,k] > H2barOut.ans;\n');
fprintf(fid, '       }\n');
fprintf(fid, '   }\n');
fprintf(fid, '}\n');

% Epsbar
fprintf(fid, 'for {i in 1..p} {\n');
fprintf(fid, '   for {j in 1..m} {\n');
fprintf(fid, '       for {k in 1..(T-1)} {\n');
fprintf(fid,...
    '            printf " %%i,%%i,%%i,%%f\\n",i,j,k,Epsbar[i,j,k] > EpsbarOut.ans;\n');
fprintf(fid, '       }\n');
fprintf(fid, '   }\n');
fprintf(fid, '}\n');

% G
fprintf(fid,' printf " %%f\\n",G >G_AzumaLP.ans;\n');

% print solver time to logfile for read-in later
fprintf(fid,...
    ' printf " %%f\\n",_solve_elapsed_time >time_AzumaLP.ans;\n');

%fprintf(fid,'\n display _total_shell_time;\n');
fprintf(fid,'\n display solve_result_num;\n');
fprintf(fid,'\n display _total_solve_time;\n');
fprintf(fid,'\n display _ampl_time;\n');

fclose(fid);

% CALL AMPL/CPLEX
system('ampl AzumaLP.run');

% read-in results to Matlab
results=readAzumaLPResults(p,m,T);

% return to original directory
cd(current)

end

function results = readAzumaLPResults(p,m,T)

% BR, 3/4/2013
% changelog: author,date,change
%{
-
%}

% NOTE - this inefficient...way to combine into single commod?

Epsbar = zeros(p,m,T-1);
fid = fopen('EpsbarOut.ans');
in = textscan(fid,'%d %d %d %f','Delimiter',',');
i=in{1};j=in{2};k=in{3};In=in{4};
numElements=length(i);
for jj = 1:numElements
    Epsbar(i(jj),j(jj),k(jj)) = In(jj);
end
fclose(fid);

H2 = zeros(m,m,T);
fid = fopen('H2Out.ans');
in = textscan(fid,'%d %d %d %f','Delimiter',',');
i=in{1};j=in{2};k=in{3};In=in{4};
numElements=length(i);
for jj = 1:numElements
    H2(i(jj),j(jj),k(jj)) = In(jj);
end
fclose(fid);

H2bar = zeros(m,m,T);
fid = fopen('H2barOut.ans');
in = textscan(fid,'%d %d %d %f','Delimiter',',');
i=in{1};j=in{2};k=in{3};In=in{4};
numElements=length(i);
for jj = 1:numElements
    H2bar(i(jj),j(jj),k(jj)) = In(jj);
end
fclose(fid);
 
G = csvread('G_AzumaLP.ans');
tSolve = csvread('time_AzumaLP.ans');

% put results into output structure
results=struct('Epsbar',Epsbar,'H2',H2,'H2bar',H2bar,'G',G,...
    'tSolve',tSolve);

end

function [p,m,T] = generateData_AzumaLP(sysC,sysParams)
% writes AzumaLP.dat for use with AMPL

% BR, 3/4/2013

% changelog: author,date,change
%{
- 3/5/2013: fixed bug with order of Phi

%}

T = sysParams.T;dt = sysParams.dt;
Pd = c2d(sysC,dt);[Ap,Bp,Cp,~] = ssdata(Pd);
% sizes:  p OUTPUTS, m INPUTS (Note this is different than automatica)
[p,~] = size(Cp);[~,m] = size(Bp);

% prep data
absCB = abs(Cp*Bp);

Phi = zeros(p,m,T);
for k = 1:T
    Phi(:,:,T-k+1) = Cp*Ap^(k-1)*Bp;
end

% NOTE: COULD MAKE MORE EFFICIENT (not listing numbers, cell array dump)
fid = fopen('AzumaLP.dat','w');

% basic dims
fprintf(fid,'param m   := %i;\n',m);
fprintf(fid,'param p   := %i;\n',p);
fprintf(fid,'param T   := %i;\n\n',T);
% control action bound
fprintf(fid,'param gamma_wv   := %f;\n\n',sysParams.gamma);

% absCB parameter
fprintf(fid,'param absCB:=\n');
for i=1:p
    for j = 1:m   
        fprintf(fid,'%i,%i    %f\n',i,j,absCB(i,j));
    end
end
fprintf(fid,';\n\n');

% Phi(i,j,k) parameter
fprintf(fid,'param Phi:=\n');
for i=1:p
    for j=1:m
        for k=1:T
            fprintf(fid,'%i,%i,%i    %f\n',i,j,k,Phi(i,j,k));
        end
    end
end
fprintf(fid,';\n\n'); 

fclose(fid);

end


function [H,Tp] = makeBlockHankel(hIn)
% makes block hankel matrix H with 3D array hIn(:,:,1:T)
% hIn(:,:,k) is an m x m matrix

[m,n,T] = size(hIn);
if(m~=n)
    disp('ERROR, nonsquare H2k')
end

Tp = floor(T/2)+1;
if(~mod(T,2))   % if T even
    hIn(:,:,T+1) = zeros(m);
    numBlocks=T+1;
else
    numBlocks=T;
end

order = (1:numBlocks)';
H = zeros(m*Tp);
for i=1:Tp
    % construct block row i
    for j=1:Tp;
        ind = order(j);
        H(m*(i-1)+1:(m*i),m*(j-1)+1:(m*j)) = hIn(:,:,ind);
    end
    order = circshift(order,-1);
end


end


% scratch attempt at linprog
%{
% construct constraint matrices for each constraint:

% constraint 1
A1 = zeros(p,numVars);
B1 = zeros(p,1);

A1(:,1) = -1*ones(p,1);
A1(:,2:end) = [ repmat(eye(p),1,(m*(T-1))), zeros(p,numVars - numEpsBar -1) ];
B1 = sum( abs(Cp*Bp), 2);

% constraint 2
A2 = zeros((m*p*(T-1))*2,numVars);

for k = 1:(T-1)

    % first prep Phi
    Phi_k = zeros(p,m*k);
    for i = 1:k
        Phi_k(:,((m*k)*(i-1)+1):(i*m*k)) = Cp*Ap^(k - i)*Bp;
    end

    % each k has (m*p*2) constraints
    % plus
    % m col strings, each with p col vecs:  m*p constraints

    epsBarCols = [repmat(zeros(m*p),1,(k-1)), kron(-eye(p),eye(m)), repmat(zeros(m*p),1,(T-1-k))];
    % each H2k has m*m vars (when expanded)
    phiHCols = [kron(Phi_k,eye(m)), repmat(zeros(m*m),1,(T-k))];
    A2thiskPlus = [zeros(m*p,1),epsBarCols, phiHCols,zeros(size(phiHCols))];
    B2thiskPlus = reshape( (Cp*Ap^k*Bp), (m*p), 1);

    % minus
    A2thiskMinus = [zeros(m*p,1),epsBarCols, -phiHCols,zeros(size(phiHCols))];
    B2thiskMinus = - B2thiskPlus;

    A2thisk = [A2thiskPlus;A2thiskMinus];
    B2thisk = [B2thiskPlus;B2thiskMinus];

end
        
A_LP = [A1;A2;A3;A4];
B_LP = [B1;B2;B3;B4];
f = [1 zeros(1,numVars-1)]; % min G

X = linprog(f,A,b);

%}