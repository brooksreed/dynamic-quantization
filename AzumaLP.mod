 # Dynamic quantizer design (Azuma 2008 TAC)
 # BR, 3/4/2013
 
 ### PARAMETERS ###
 param m;
 param p;
 param T;
 param gamma_wv;
 param absCB{1..p,1..m};
 param Phi{1..p,1..(m*(T-1)),1..T};
 
 ### VARIABLES ###
 var G;
 var H2{1..m,1..m,1..T};
 var H2bar{1..m,1..m,1..T} >=0;
 var Epsbar{1..p,1..m,1..(T-1)} >=0; 
 var innerII{1..p,1..m,1..(T-1)};  # auxiliary var for constraint II
 
 ### OBJECTIVE ###
 minimize Cost:	G;
 
 ### CONSTRAINTS ###
 
 subject to I{i in 1..p}:
	sum{j in 1..m} absCB[i,j] 
	+ sum{j in 1..m} sum{k in 1..(T-1)} Epsbar[i,j,k] <= G;
	
 subject to insideII{i in 1..p,j in 1..m,k in 1..(T-1)}:
	innerII[i,j,k] = sum{kk in 1..(k)} sum{mm in 1..m} Phi[i,mm,(T-kk+1)]*H2[mm,j,(k-kk+1)];
	
	#Phi = zeros(p,m,T);for k = 1:T; Phi(:,:,T-k+1) = Cp*Ap^(k-1)*Bp; end
	
 subject to II_p{i in 1..p,j in 1..m,k in 1..(T-1)}:
	Phi[i,j,(T-k)] + innerII[i,j,k] <= Epsbar[i,j,k];
	
 subject to II_m{i in 1..p,j in 1..m,k in 1..(T-1)}:
	Phi[i,j,(T-k)] + innerII[i,j,k] >= -Epsbar[i,j,k];
 
 subject to III{i in 1..m}:
	1 + sum{k in 1..(T)} sum{j in 1..m} H2bar[i,j,k] <= gamma_wv;
	
 subject to IV_p{i in 1..m,j in 1..m, k in 1..(T)}:
	-H2bar[i,j,k] <= H2[i,j,k];

 subject to IV_m{i in 1..m,j in 1..m, k in 1..(T)}:
	H2bar[i,j,k] >= H2[i,j,k];
	
 