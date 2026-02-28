clear
tic
global sdownc supc annic creac nn gamma kappa gammap H sz
%% parameter setting and basic operator definition
w0 = 1; wc = 1; g = 0.01; nn = 10; beta = 1; gamma = 1e-4; kappa = 1e-4; gammap = 1e-4;
NN = 50;% number of rounds
n = sqrt(linspace(0,nn-1,nn));
anni = full(spdiags(n',1,nn,nn));
crea = anni';
annic = kron(eye(2),anni);
creac = kron(eye(2),crea);
sup = [0 1;0 0];
sdown = sup';
supc = kron(sup,eye(nn));
sdownc = kron(sdown,eye(nn));
sz = kron([1,0;0,-1],eye(nn));
Mg = [0,0;0,1];
Me = [1,0;0,0];
%% Hamiltonian setting
H0=kron([1,0;0,0],eye(nn))+wc*kron(eye(2),crea*anni);
Hi=g*kron([0,1;0,0],anni);
HI=Hi+Hi';
H=0*H0+HI;
%% initial state setting
inis = [1;0];% initial state of the qubit
tn=1;
aa = sqrt(tn);% generating single fock state
% aa=factorial(tn)^(1/(2*tn));% generating Fock sate superposition
inic = zeros(nn,1);
% psi0=kron(inis,inic);
for kk = 1:length(inic)
    inic(kk) = exp(-0.5*aa^2)*aa^(kk-1)/sqrt(factorial(kk-1));% coherent state
end

% Thermal state
% inic=zeros(nn,1);
% rho0ctemp=diag(exp(-beta*wc*linspace(0,nn-1,nn))); % thermal state
% rho0c=rho0ctemp/trace(rho0ctemp);

% initial full density matrix
rho0c = inic*inic';
rho0 = kron(inis*inis',rho0c);
rhotf = zeros(nn,nn);
%% target state
% single Fock state
target = zeros(nn,1);
target(tn+1) = 1;
rhotarget = target*target';

% superposed Fock
% target = zeros(nn,1);
% target(1) = 1;
% target(tn+1) = 1;
% rhotargetp = (target*target')/trace(target*target');
% target(1) = 1;
% target(tn+1) = -1;
% rhotargetm = (target*target')/trace(target*target');
% rhotarget = rhotargetm+rhotargetp;
% rhotargetm2 = rhotargetm;
% rhotargetp2 = rhotargetp;
%%
Pn = zeros(1,nn);
Pg=zeros(1,NN+1);
Pg(1) = 1;
Pt = zeros(1,NN);

Omegatn = sqrt(g^2*(tn+1));
tau = 1*pi/Omegatn;
rhot = rho0;
 for ii = 1:NN
        yt = reshape(rhot,[4*nn*nn,1]);
        [tt,yy] = ode45(@Lindbladcav,linspace(0,tau,1000),yt);
        rhot = reshape(yy(length(tt),:).',[2*nn,2*nn]);
        rhottemp = kron(Me,eye(nn))*rhot*kron(Me,eye(nn));
        rhot = rhottemp/trace(rhottemp);
        Pg(ii+1) = trace(rhottemp)*Pg(ii);
        Pt(ii) = trace(rhot*kron([1 0;0 0],rhotarget));
        % PP1(jj,ii)=trace(rhot*kron([1 0;0 0],rhotargetp1))/trace(rhot*kron([1 0;0 0],rhotargett));
        % PM1(jj,ii)=trace(rhot*kron([1 0;0 0],rhotargetm1))/trace(rhot*kron([1 0;0 0],rhotargett));
        % disp(Pt(jj,ii));
 end
 for kk = 1:nn
        tar = zeros(nn,1);
        tar(kk) = 1;
        Pn(kk) = trace(rhot*kron(eye(2),tar*tar'));
 end
 toc
%%
figure(1),
plot(1:1:NN,Pt), hold on
%%
figure(2),
bar(0:1:nn-1,Pn(1,:))
