%% MAE509 LMI Wikibook Contribution
% Josh Stoffel
clear all; clc
opts = sdpsettings('solver','mosek','verbose',0);
tol = 1e-6;

%% F-16 Linearized Longitudinal Model

A = [-1.9311e-2,    8.8157, -3.217e1, -5.7499e-1; 
     -2.5389e-4,   -1.0189,        0,  9.0506e-1;
              0,         0,        0,          1;
     2.9465e-12, 8.2225e-1,        0,    -1.0774];
 
B1 = [ 1.737e-1;
    -2.1499e-3;
             0;
    -1.7555e-1];

B2 = [ 1.737e-1;
    -2.1499e-3;
             0;
    -1.7555e-1];

C1 = [0, 5.729578e1, 0,          0;
     0,          0, 0, 5.729578e1];

C2 = [0.003981, 15.88, 0,       1.481;
     0,          0, 0, 0];
 
D11 = [0; 0];

D12 = [0; 0];

D21 = [0; 1];

D22 = [0.03333; -1];


%% H-infinity Dynamic Output Feedback Controller Synthesis
tol = 0.000001;   
ns = size(A,1);   % number of states
nc = size(B2,2);  % number of actuators
nm = size(C2,1);  % number of sensors
nd = size(B1,2);  % number of external inputs
no = size(C1,1);  % number of regulated outputs

% Optimization Variables
gamma = sdpvar;         
X1 = sdpvar(ns);
Y1 = sdpvar(ns);
An = sdpvar(ns,ns,'full');
Cn = sdpvar(nc,ns,'full');
Dn = sdpvar(nc,nm,'full');
Bn = sdpvar(ns,nm,'full');

% Linear Matrix Inequality Constraints
cons = [X1 >= tol*eye(ns), Y1 >= tol*eye(ns)];
cons = [cons, [X1 eye(ns); eye(ns) Y1] >=0 ];
LMI = [A*Y1+Y1*A'+B2*Cn+Cn'*B2'  (A'+An+(B2*Dn*C2)')'        B1+B2*Dn*D21           (C1*Y1+D12*Cn)'; 
       A'+An+(B2*Dn*C2)'         X1*A+A'*X1+Bn*C2+C2'*Bn'    X1*B1+Bn*D21           (C1+D12*Dn*C2)'  ;
       (B1+B2*Dn*D21)'           (X1*B1+Bn*D21)'             -gamma*eye(nd)          (D11+D12*Dn*D21)'  ;
       C1*Y1+D12*Cn              C1+D12*Dn*C2                D11+D12*Dn*D21         -gamma*eye(no)];

cons = [cons, LMI<=0];

% Minimize Hinf Norm
optimize(cons,gamma,opts)
gamman = value(gamma);
X1n = value(X1);
Y1n = value(Y1);
Ann = value(An);
Bnn = value(Bn); 
Cnn = value(Cn);
Dnn = value(Dn);
K1=[Ann Bnn; Cnn Dnn]-[X1n*A*Y1n zeros(ns,nm); zeros(nc,ns) zeros(nc,nm)];

% Reverse variable substitution
Y2n = eye(ns);
X2n = inv(Y2n)*(eye(ns)-X1n*Y1n);
K2 = inv([X2n X1n*B2;zeros(nc,ns) eye(nc)])*K1*inv([Y2n' zeros(ns,nm); C2*Y1n eye(nm)]);
Ak2 = K2(1:ns,1:ns);Bk2=K2(1:ns,(ns+1):(ns+nm));Ck2=K2((ns+1):(ns+nc), 1:ns);Dk2=K2((ns+1):(ns+nc), (ns+1):(ns+nm));
Dk = inv(eye(nc)+Dk2*D22)*Dk2;
Bk = Bk2*(eye(nm)-D22*Dk);
Ck = (eye(nc)-Dk*D22)*Ck2;
Ak = Ak2-Bk*inv(eye(nm)-D22*Dk)*D22*Ck;

% Lower LFT Closed Loop System
plant = ss(A,[B1 B2],[C1;C2],[D11 D12; D21 D22]);
controller = ss(Ak,Bk,Ck,Dk);
sys_cl = lft(plant,controller);
norm(sys_cl,Inf);