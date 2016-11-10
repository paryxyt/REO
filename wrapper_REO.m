 %% test example
 clear
 clc
 input_data=@input_data_RE; % specify the input data here
 [dims,A,b,c,Kparams]=input_data(); 
 [y, nu, obj] = call_solver_REO(dims,A,b,c,Kparams)


 %% signomial
 % wrapper for relative entropy
 A=[0 10.2000         0         0    1.5089    1.0857    1.0459  
    0      0    9.8000         0    1.0981    1.9069    0.0492 
   0      0         0    8.2000    1.3419    1.6192    1.6245 ];
 
     
 c=[0 10.0000   10.0000   10.0000  -14.6794   -7.8601    8.7838]';
 %f_sage=-.9747;
 
 %A=[ 0   10.2000         0         0    1.3071   -3.4526    0.2864
 %        0         0    9.8000         0   -2.2097   -1.6080   -0.6931
 %        0         0         0    8.2000   -1.3975    2.9699   -0.7666];
%c= [   0   10.0000   10.0000   10.0000   -1.0092   15.2330   22.8559  ]';
 
 % variable of REO problem are [v; tau_1, ..., tau_n]
 % each tau_i is m dimensional
 n=7; %number of terms;
 m=3; %number of variables
 
 sz=nchoosek(n,2);
 B=zeros(n+m*n,2*sz);
 C=zeros(n+m*n,2*sz);
 D=zeros(n+m*n,2*sz);
 
 idx=0;
 for i=1:n
     for j=1:n
         if i~=j
             idx=idx+1;
             B(i,idx)=1;
             C(j,idx)=1;
             D(n+(i-1)*m+1:n+(i-1)*m+m,idx)=-(A(:,i)-A(:,j));
         end
     end
 end
 
 % create dims, A, b , c0, Kparams
 dims=[n+m*n 1 2*sz]; %number of variables in REO problem, number of equality constraints,  number RE constraints 
  
 c0=zeros(n+m*n,1);
 c0(1:n)=c;
 
 A=zeros(1,n+m*n);
 A(1)=1;
b=1;
 
Kparams=[B C D];

[y, nu, obj] = call_solver_REO(dims,A,b,c0,Kparams)

 