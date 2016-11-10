function [dims,A,b,c,Kparams]=input_data3()
% Minimize:    c'*y
% subject to:  Bi'*y log(Bi'*y/Ci'*y)+ Di'*y <= 0    (*)


% defines problem parameters for RE cone
n1=2; % size of y
m1=1; % number of RE constraints
m=0; %number of equality constraints
dims=[n1 m m1];

% there are m1 inequalities of the form (*)
B=zeros(n1,m1);
C=zeros(n1,m1);
D=zeros(n1,m1);

% B=[ 1 0 0 0 0 0; 
%     0 1 0 0 0 0;
%     0 0 1 0 0 0];
B=[1 0];
B=B';
C=[0 1];
% C=[1 .1 .2 0 0 0;
%   .1  1  0 0 0 0;
%   .2  0 1 0 0 0 ];
C=C';
%D=-rand(1,2);
D=[0 -1];
% D=[0 0 0   -.1    0   0;
%    0 0 0     0  -.1   0;
%    0 0 0     0    0 -.1];
D=D';

%c=rand(n1,1);
c=[1; 1];

Kparams=[B C D];

A=[];
b=[];