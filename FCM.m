function [Iu,U,P,norm_P]=FCM(Data,P0,M,epsm)
% Fuzzy c-means clustering (FCM) : began to iterate
% from the initial clustering center
% [U,P,Dist,Cluster_Res,Obj_Fcn,iter] = fuzzycm2(Data,P0,plotflag,M,epsm)
% inout: Data,plotflag,M,epsm: see fuzzycm.m
% P0:initial clustering center
% output: U,P,Dist,Cluster_Res,Obj_Fcn,iter: see fuzzycm.m
% See also: fuzzycm
[N,S] = size(Data); m = 2/(M-1); iter = 0;
C=size(P0,1);Dist(C,N)=0;U(C,N)=0;P(C,S)=0;
% The iteration algorithm of FCM
while true
    % Iterative counter
    iter=iter+1;
%     fprintf('*');
    % Calculate or update divided matrix U
   
    for i=1:C
        Dist(i,:)= sqrt(sum((Data-ones(N,1)*P0(i,:)).^2,2))';
    end
    U=1./(Dist.^m.*(ones(C,1)*sum(Dist.^(-m))));
    U(Dist==0)=1;

    % update clustering center p
    Um=U.^M;
    P=Um*Data./(ones(S,1)*sum(Um'))';

    % iterative stop condition of FCM algorithm
    if norm(P-P0,Inf)<epsm
       %fprintf('Obj_Fcn=%d,', Obj_Fcn(iter));
       break
    end
    
    if iter>100
%         iter>500
        
        break
    end
    P0=P;
end
norm_P=norm(P-P0,Inf);
fprintf('FCM\n');
[o,Iu]=max(U);