%% needs symbolic toolbox

syms x1 x2 x3 s
f=-log(-x2*log(x2/x3)-x1+s)-log(x2)-log(x3);
df_dx1=diff(f,x1)
df_dx2=diff(f,x2)
df_dx3=diff(f,x3)
df_ds=diff(f,s)

dx1x1=diff(df_dx1,x1)
dx1x2=diff(df_dx1,x2)
dx1x3=diff(df_dx1,x3)
dx2x2=diff(df_dx2,x2)
dx2x3=diff(df_dx2,x3)
dx3x3=diff(df_dx3,x3)
dx1s=diff(df_dx1,s)
dx2s=diff(df_dx2,s)
dx3s=diff(df_dx3,s)
dss=diff(df_ds,s)