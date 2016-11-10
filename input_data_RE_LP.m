function [dims,A,b,c,Kparams]=input_data_RE_LP()
rng(2011,'twister')

% Minimize:    c'*y
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= 0    (*)


% defines problem parameters for RE cone
n1=6; % size of y
n2=2; % size of z
m=4; %number of equality constraints
m1=3; % number of RE constraints
dims=[n1 n2 m m1];

% there are m1 inequalities of the form (*)
B=zeros(n1,m1);
C=zeros(n1,m1);
D=zeros(n1,m1);

 B=[ 1 0 0 0 0 0; 
     0 1 0 0 0 0;
     0 0 1 0 0 0];
B=B';
C=[1 .1 .2 0 0 0;
   .1  1  0 0 0 0;
   .2  0 1 0 0 0 ];
C=C';
D=[0 0 0   -.1    0   0;
    0 0 0     0  -.1   0;
    0 0 0     0    0 -.1];
D=D';


Kparams=[B C D];

A=rand(m,n1+n2);
b=rand(m,1);
c=rand(n1+n2,1);
