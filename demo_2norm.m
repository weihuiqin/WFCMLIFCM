clear
clc
close all
tic
%load('ORL_mtv.mat');
load('MF_mtv.mat');
cluster = length(unique(gt));
threshold = 1e-8;
delta=1000000;%??????????????
zeta=1000;
ss=6;
N=length(gt);
P=length(X);
U=cell(P,1);
center=cell(P,1);
center1=cell(P,1);
S=cell(P,1);
for p=1:P
   S{p}=size(X{p},1); 
end
m_degree=[1.2,1.2,1.2,1.2,1.2,1.2];
for p=1:P
   X{p} = X{p}./(repmat(sqrt(sum(X{p}.^2,1)),size(X{p},1),1)+10e-10);
end

parfor p=1:P
[U{p,1},nmi1{p,1},nmi2{p,1},acc1{p,1},acc2{p,1},center1{p,1}] = k_fcm(X{p},gt,cluster,m_degree(p),threshold,ss);
end
load('UUU.mat');

for p=1:P
U{p,1}=UU;
end

for p=1:P
center{p,1}=center1{p,1};
end

Dist1=cell(P,1);
Dist_star(cluster,N)=0;
for p=1:P
   alpha(p)=1/P; %??????weight
end
gamma=0.1;
iter=0;
Isconverg = 0;
m_degree=2;
m=1/(m_degree-1);
M=2;
epsm=threshold;
[U1,U_star] = norm2(X, cluster,1,M,epsm,center,delta,zeta,U,N,P,S);
U=U1;
[o,Iu]=max(U_star);
[nmi,acc]=measure(Iu',gt);
U_tensor = cat(3, U{:,:});
Tensor_U=tensor(U_tensor);
Error_Threshold=1e-8;
Max_iterations=100;
m_degree=1.2;
  for i=cluster
         for j=cluster
             for k=P
             Rank_A=[i,j,k];
             [Core_Tensor_A,Singular_Factors_A,Singular_FA]=Decompose_Tensor_Coupled_HOSVD_iteratively(Tensor_U,Rank_A,Error_Threshold,Max_iterations);
             A1_n=Singular_Factors_A{1,1};
             % data=bsxfun(@rdivide,bsxfun(@minus,meas,min(meas)),max(meas)-min(meas));
             [idx3,P1] = kmeans1(A1_n,cluster);
             [nmi3(i,j,k),acc3(i,j,k)]=measure(idx3,gt);
  
             [idx4,U1]=FCM(A1_n,P1,m_degree,threshold);
             idx4 = idx4'; 
             [nmi4(i,j,k),acc4(i,j,k)] = measure(idx4,gt);
  
              A1_nn=Singular_FA{1,1};
              [idx5,P2] = kmeans1(A1_nn,cluster);
              [nmi5(i,j,k),acc5(i,j,k)]=measure(idx5,gt);

              [idx6,~]=FCM(A1_nn,P2,m_degree,threshold);
              idx6 = idx6'; 
              [nmi6(i,j,k),acc6(i,j,k)] = measure(idx6,gt);
              
              A1_nn=Singular_FA{1,1};
              [idx7,P3] = kmeans1(UU,cluster);
              [nmi7(i,j,k),acc7(i,j,k)]=measure(idx7,gt);

              [idx8,~]=FCM(UU,P3,m_degree,threshold);
              idx8 = idx8'; 
              [nmi8(i,j,k),acc8(i,j,k)] = measure(idx8,gt);
  
             end
         end
  end
   a1=max(max(max(acc3)));
   a2=max(max(max(acc4)));
   a3=max(max(max(acc5)));
   a4=max(max(max(acc6)));
   a5=max(max(max(nmi3)));
   a6=max(max(max(nmi4)));
   a7=max(max(max(nmi5)));
   a8=max(max(max(nmi6)));
   a9=max(max(max(acc7)));
   a10=max(max(max(acc8)));
   a11=max(max(max(nmi7)));
   a12=max(max(max(nmi8)));
   toc