function [U1,U_star] = norm2(X, cluster,plotflag,M,threshold,center,delta,zeta,U,N,P,S)
%N=length(gt);
%P=length(X);
Dist1=cell(P,1);
Dist_star(cluster,N)=0;
for p=1:P
   alpha(p)=1/P; 
end
gamma=0.1;
center1=cell(P,1);
% eta = 10e-5;
% eta = 10e-9;
% eta = 10e-7; max_eta = 10e10; pho_eta = 2;
%while (Isconverg == 0)
if nargin<5
    threshold=1.0e-3;
end
if nargin<4
    M=2;
end
if nargin<3
    plotflag=0;
end
max_iter=1;
min_impro=threshold;
expo=M;
m_degree=2;
m=1/(m_degree-1);
for iter=1:max_iter
 % Iterative counter
 %iter=iter+1;
 fprintf('----processing iter %d--------\n', iter+1);
 % update U^p
 %%%%%%%%%%% calculate local membership values%%%%%%%%%%%%%%%%%%
   for p=1:P
    for i=1:cluster
        Dist(i,:)= sqrt(sum((X{p}'-ones(N,1)*center{p,1}(i,:)).^2,2))';
    end
    R=2*(Dist.^2)+2*delta*(P-1)^2/(P^2);
    %R=2*(Dist.^2)+eta*T*ones(cluster,N);
    R1=R.^(-1);
    R1=R1';
    S1=0;
    for pp=1:P
    S1=S1+U{pp};
    end
    S2=S1-U{p};
    S2=2*delta*(P-1)*S2/(P^2);
    tmp=(sum(S2.*R1,2)-1)./(sum(R1,2));
    U{p,1}=(-tmp*ones(1,cluster)+S2).*R1;
    U{p,1}(R==0)=1;
   end
   %%%%%%%%%%% calculate global membership values without weights%%%%%%%%%%%%%%%%%%
for p=1:P
  for i=1:cluster
        Dist1{p,1}(i,:)= (sum((X{p}'-ones(N,1)*center{p,1}(i,:)).^2,2))';
  end
     Dist_star=Dist_star+(alpha(p)^gamma)*Dist1{p,1};
end
   U_star=1./(Dist_star.^m.*(ones(cluster,1)*sum(Dist_star.^(-m))));
   U_star(Dist_star==0)=1;
   % update center^p
   %%%%%%%%%%% calculate the centers %%%%%%%%%%%%%%%%%%
   for p=1:P
    Um3=U{p,1}.^expo+zeta*(U_star'.^expo);
    center1{p,1}=Um3'*X{p}'./(ones(S{p},1)*sum(Um3))';
   end
   %%%%%%%%%%%%%%coverge condition %%%%%%%%%%%
   %Isconverg = 1; 
   
%    for p=1:P
%        if (norm(center{p,1}-center1{p,1},Inf)>threshold)
%           history.norm_center=norm(center{p,1}-center1{p,1},Inf);
%           %Isconverg = 0;
%        end
%       center{p,1}=center1{p,1};
%    end
if nargout>4 | plotflag
        fprintf('*');
    end  
%    if iter > 1,
%        for p=1:P
% 		if norm(center{p,1} - center1{p,1},Inf) < min_impro, break; end,
%        end
%    end
%     for p=1:P
%    center{p,1}=center1{p,1};
%     end
%     if iter > 1,
%        for p=1:P
% 		if norm(center{p,1} - center1{p,1},Inf) < min_impro, break; end,
%         norm(center{p,1} - center1{p,1},Inf)
%        end
%     end
   if iter > 1,
       %for p=1:P
		if norm(center{6,1} - center1{6,1},Inf) < 1.0e-4, break; end,
        norm(center{6,1} - center1{6,1},Inf)
       %end
   end
    for p=1:P
   center{p,1}=center1{p,1};
    end
%    if (iter>2)
%         Isconverg  = 1;
%    end
   
end

for p=1:P
U{p,1}=abs(U{p,1});
end
U1=U;

