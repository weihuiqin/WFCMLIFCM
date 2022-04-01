function [Menbership,nmi1,nmi2,acc1,acc2,P] = k_fcm(X, gt,cluster,m_degree,threshold,ss)
rand('state',ss-1);
 for iter =1:1;
%%----------------------------------kmeans--------------------------
      [idx, P] = kmeans1(X',cluster);
      [nmi1(1,iter),acc1(1,iter)]= measure(idx,gt);
%     ari1(1,iter)=adjrand(dx_REG,ground_truth);

   %%----------------------------------fcm--------------------------
      [idx,U]=FCM(X',P,m_degree,threshold);
%     [idx,U,Obj_fcn(iter)]=FCM1(X{k}',P,1,m_degree,threshold);
%     [~,idx1]=min(Obj_fcn);
      A_t(:,:,iter) = U';
      idx = idx'; 
      [nmi2(1,iter),acc2(1,iter)] = measure(idx,gt);
%     ari2(1,iter)=adjrand(dx_REG,ground_truth);    
 end
      Menbership=A_t;
end
    