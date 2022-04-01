clear
clc
close all

load('ORL_mtv.mat');
cluster = length(unique(gt));
% m_degree = 1.01;% 预设值需要考虑
threshold = 1e-8;
% delta=1/3;%预设值需要考虑
delta=0.0000001;%预设值需要考虑
ss=1;
T=3;
N=length(gt);
P=length(X);
U=cell(P,1);
sX = [N, cluster, P];
nmi1=cell(P,1);
nmi2=cell(P,1);
acc1=cell(P,1);
acc2=cell(P,1);
center=cell(P,1);
S=cell(P,1);
m_degree=[2,2,2,2,2,2];
for p=1:P
   X{p} = X{p}./(repmat(sqrt(sum(X{p}.^2,1)),size(X{p},1),1)+10e-10);
end

for p=1:P
W{p}=zeros(N,cluster);
S1{p}=zeros(N,cluster);
S{p}=size(X{p},1);
end
S1_tensor = cat(3, S1{:,:});

parfor p=1:P
[U{p,1},nmi1{p,1},nmi2{p,1},acc1{p,1},acc2{p,1},center1{p,1}] = k_fcm(X{p},gt,cluster,m_degree(p),threshold,ss);
end

m_degree = 2;%预设值需要考虑
Isconverg = 0;
ModCount = 3;%tensor 展开时3种情况
for v=1:ModCount
    para_ten{v} = delta;
end
for i=1:ModCount
    H_tensor{i} = [];%复制3份W_tensor到WT 
end
 W_tensor = cat(3, W{:,:});
 
for i=1:ModCount
    WT{i} = W_tensor;%复制3份W_tensor到WT 
end

iter=0;
% eta = 10e-5;
% eta = 10e-9;
eta = 10e-7; max_eta = 10e10; pho_eta = 2;
while(Isconverg == 0)
 fprintf('----processing iter %d--------\n', iter+1);
    % update H^t
   U_tensor = cat(3, U{:,:});
   for umod=1:ModCount
        H_tensor{umod} = updateH_tensor(WT{umod},U,sX,eta,para_ten,P,umod);%需要多次检查
        WT{umod} = WT{umod}+eta*(U_tensor-H_tensor{umod});
   end
   
   for p=1:P
    % update center^p
    Um=U{p,1}.^m_degree;
    center{p,1}=Um'*X{p}'./(ones(S{p},1)*sum(Um))';
    % update U^p
    for i=1:cluster
        Dist(i,:)= sqrt(sum((X{p}'-ones(N,1)*center{p,1}(i,:)).^2,2))';%需要检查
    end
    R=2*(Dist.^2)+eta*T*ones(cluster,N);%需要检查
    R1=R.^(-1);
    R1=R1';
    
    for  umod=1:ModCount
       S1_tensor=S1_tensor+abs(H_tensor{umod}-WT{umod});%%%%%%%%%需要检查  强行设置成为正数
    end
    S1_tensor=eta*S1_tensor;
    S2=S1_tensor(:,:,p);%%%%%%%%%需要检查
    tmp=(sum(S2.*R1,2)-1)./(sum(R1,2));%%%重点检查
    U{p,1}=(-tmp*ones(1,cluster)+S2).*R1;%重点检查
    U{p,1}(R==0)=1;
   end
 %% coverge condition
    Isconverg = 1;    
   for p=1:P
       if (norm(center{p,1}-center1{p,1},Inf)>threshold)
          history.norm_center=norm(center{p,1}-center1{p,1},Inf);
          Isconverg = 0;
       end
      center1{p,1}=center{p,1};
   end
   for umod=1:ModCount
       for p=1:P
       if (norm(U_tensor(:,:,p)-H_tensor{umod}(:,:,p),inf)>threshold)
           history.norm_U_H=norm(U_tensor(:,:,p)-H_tensor{umod}(:,:,p),inf);
           Isconverg = 0;
       end
       end
   end
 if (iter>2)
        Isconverg  = 1;
 end
    iter = iter + 1;
    eta = min(eta*pho_eta, max_eta);
end
for p=1:P
U{p,1}=abs(U{p,1});
end
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
  
%               ar=1.1;
%               [idx7,U2]=EAWFCM(A1_n,cluster,m_degree,ar);
%               idx7 = idx7'; 
%               [nmi7(i,j,k),acc7(i,j,k)] = measure(idx7,gt);
  
  
              A1_nn=Singular_FA{1,1};
              [idx5,P2] = kmeans1(A1_nn,cluster);
              [nmi5(i,j,k),acc5(i,j,k)]=measure(idx5,gt);

              [idx6,~]=FCM(A1_nn,P2,m_degree,threshold);
              idx6 = idx6'; 
              [nmi6(i,j,k),acc6(i,j,k)] = measure(idx6,gt);
  
%               m_degree=1.2356789;
%               ar=9;
%               [idx8,U3]=EAWFCM(A1_nn,cluster,m_degree,ar);
%               idx8 = idx8'; 
%               [nmi8(i,j,k),acc8(i,j,k)] = measure(idx8,gt);
             end
         end
  end
   max(max(max(acc3)))
   max(max(max(acc4)))
   max(max(max(acc5)))
   max(max(max(acc6)))
   max(max(max(nmi3)))
   max(max(max(nmi4)))
   max(max(max(nmi5)))
   max(max(max(nmi6)))