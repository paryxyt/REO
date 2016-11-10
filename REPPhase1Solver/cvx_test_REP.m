clc
rng(2010,'twister')

%Solver for REPs using Primal only method
% Minimize:    sum(s1)
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= s1_i    (*)
%              A1*y + A2*z = b                         (**)
%              z >= 0
% Calls grad_f, PrimalOnlyNewtonResidual, infeasibleNewton



%% test example data 1 (uncomment)
% n1=6; % size of y
% n2=2; % size of z
% m1=3; % size of s1
% m2=2; % num equalities of form (**)
% dims=[n1 n2 m1 m2];
% 
% 
% % there are m1 inequalities of the form (*)
% B=zeros(n1,m1);
% C=zeros(n1,m1);
% D=zeros(n1,m1);
% 
% B=[ 1 7 4 0 0 0; 
%     3 1 6 0 0 0;
%     6 1 8 0 0 0];
% B=B';
% 
% C=[1 3 4 0 0 0;
%     2 8  1 0 0 0;
%     9 1 2 0 0 0 ];
% C=C';
% 
% D=[0 0 0   -.1    0   0;
%    0 0 0     0  -.1   0;
%    0 0 0     0    0 -.1];
% D=D';
% 
% A1=rand(m2,n1); % there are m2 equalities of the form (**)
% A2=rand(m2,n2); % the dimension of z is n2
% b=zeros(m2,1);
% %b=rand(m2,1);
% 
% A=[A1 A2 zeros(m2,m1)];
% c=zeros(n1+n2+m1,1);
% c(n1+n2+1:n1+n2+m1,1)=1;


%% test example data 2 (simple example)
n1=3; % size of y
n2=1; % size of z
m1=1; % size of s1
m2=1; % num equalities of form (**)
dims=[n1 n2 m1 m2];


% there are m1 inequalities of the form (*)
B=zeros(n1,m1);
C=zeros(n1,m1);
D=zeros(n1,m1);

B=[1 0 0];
C=[0 1 0];
D=[0 0 1];

A1=[1 1 1];
A2=1;
b=0;
% B=[ 1 7 4 0 0 0; 
%     3 1 6 0 0 0;
%     6 1 8 0 0 0];
 B=B';
% 
% C=[1 3 4 0 0 0;
%     2 8  1 0 0 0;
%     9 1 2 0 0 0 ];
 C=C';
% 
% D=[0 0 0   -.1    0   0;
%    0 0 0     0  -.1   0;
%    0 0 0     0    0 -.1];
 D=D';

% A1=rand(m2,n1); % there are m2 equalities of the form (**)
% A2=rand(m2,n2); % the dimension of z is n2
% b=zeros(m2,1);
%b=rand(m2,1);

A=[A1 A2 zeros(m2,m1)];
c=zeros(n1+n2+m1,1);
c(n1+n2+1:n1+n2+m1,1)=1;

cvx_begin
variable y(n1,1)
variable z(n2,1)
variable s1(m1,1)
%variable s2(n2,1)
X=[y;z;s1];

minimize sum(s1)
subject to
rel_entr(B(:,1)'*y,C(:,1)'*y) + D(:,1)'*y <= s1(1)
%rel_entr(B(:,2)'*y,C(:,2)'*y) <= D(:,2)'*y + s1(2)
%rel_entr(B(:,3)'*y,C(:,3)'*y) <= D(:,3)'*y + s1(3)
A*X==b
z>=0
s1>=-1
cvx_end

%check order of inequalities