rng(2011,'twister')
clear
%% LP example
input_data=@input_data_LP; % specify the input data here
[dims,A,b,c,Kparams]=input_data();
n=dims(1);
m=dims(2);
 
 cvx_begin
 variable x(n,1)
 minimize transpose(c)*x
 subject to 
 A*x==b
 x>=0
 cvx_end

%% SDP example
clear
rng(2011,'twister')

input_data=@input_data_SDP; % specify the input data here
[dims,A,b,c,Kparams]=input_data();
n=dims(1);
m=dims(2);
p=dims(3);


 cvx_begin sdp
 variable X(p,p) symmetric
 minimize transpose(c)*vec(X)
 subject to 
 A*vec(X)==b
 X>=0
 cvx_end


%% RER example
input_data=@input_data4; % specify the input data here
[dims,A,b,c,K_params]=input_data();
n=dims(1);
m=dims(2);
m1=dims(3);

B=K_params(:,1:m1);
C=K_params(:,m1+1:2*m1);
D=K_params(:,2*m1+1:3*m1);

cvx_begin
variable y(n,1);
minimize transpose(c)*y
subject to
%rel_entr(B(:,1)'*y,C(:,1)'*y) + D(:,1)'*y <= 0
rel_entr(B'*y,C'*y) + D'*y <= 0
A*y==b
cvx_end

c
A
y

%% RER example Phase 1
rng(2011,'twister')
%input_data=@input_data_RE_Phase1; % specify the input data here; % specify the input data here
input_data=@input_data_RE;
[dims,A,b,c,K_params]=input_data();
n=dims(1);
m=dims(2);
m1=dims(3);

B=K_params(:,1:m1);
C=K_params(:,m1+1:2*m1);
D=K_params(:,2*m1+1:3*m1);

cvx_begin
variable y(n,1);
minimize transpose(c)*y
subject to
%rel_entr(B(:,1)'*y,C(:,1)'*y) + D(:,1)'*y <= 0
rel_entr(B'*y,C'*y) + D'*y <= 0
A*y==b
cvx_end

c
A
y


%% RER + LP example Phase 1
clear
clc
rng(2010,'twister')
%input_data=@input_data_RE_Phase1; % specify the input data here; % specify the input data here
input_data=@input_data_RE_LP; % specify the input data here
[dims,A,b,c,Kparams]=input_data();
 %n=dims(1)+dims(4)+dims(2);
 n=dims(1)+dims(2);
 n1=dims(1); % number of RE variables
 n2=dims(2); % LP variables
 m=dims(3); % number of equality constraints
 m1=dims(4);  % number of RE constraints
 dims_phase1=[n m m1 n1 n2];

B=Kparams(:,1:m1);
C=Kparams(:,m1+1:2*m1);
D=Kparams(:,2*m1+1:3*m1);

A_phase1 = [A(:,1:n1) zeros(m,m1) A(:,n1+1:n1+n2)];
 c_phase1=zeros(n1+n2+m1,1);
 c_phase1(n1+1:n1+m1)=1; 

A=A_phase1;
%c=c_phase1;
cvx_begin
variable y(n,1);
y1=y(1:n1);
s=y(n1+1:n1+m1);
y2=y(n1+m1+1:n1+m1+n2);

minimize sum(s)
subject to
%rel_entr(B(:,1)'*y,C(:,1)'*y) + D(:,1)'*y <= 0
rel_entr(B'*y1,C'*y1) + D'*y1 <= 0
A*y==b
y2>=0
cvx_end

c
A
y

%% RER + LP example Phase 2
clear
clc
rng(2010,'twister')
%input_data=@input_data_RE_Phase1; % specify the input data here; % specify the input data here
input_data=@input_data_RE_LP; % specify the input data here
[dims,A,b,c,Kparams]=input_data();
 n=dims(1)+dims(2);
 n1=dims(1); % number of RE variables
 n2=dims(2); % LP variables
 m=dims(3); % number of equality constraints
 m1=dims(4);  % number of RE constraints
 dims_phase1=[n m m1 n1 n2];

B=Kparams(:,1:m1);
C=Kparams(:,m1+1:2*m1);
D=Kparams(:,2*m1+1:3*m1);

cvx_begin
variable y(n,1);
y1=y(1:n1);
y2=y(n1+1:n1+n2);

minimize transpose(c)*y
subject to
%rel_entr(B(:,1)'*y,C(:,1)'*y) + D(:,1)'*y <= 0
rel_entr(B'*y1,C'*y1) + D'*y1 <= 0
A*y==b
y2>=0
cvx_end

c
A
y