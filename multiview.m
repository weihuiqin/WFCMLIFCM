function [A_t] = multiview(X, gt,cluster,m_degree,threshold,ss)
V = length(X);
Menbership=cell(V,1);
 parfor k=1:V
   Menbership{k,1} = k_fcm(X{k},gt,cluster,m_degree,threshold,ss);
 end
 for k=1:V
 A_t(:,:,(1+(k-1)*10):10*k)=Menbership{k,1};
 end
end

% final_acc11=sum(acc_11)/iter;
% final_nmi11=sum(nmi_11)/iter;
% fprintf('kmeans_acc1\n%.4f\n',final_acc11);
% fprintf('kmeans_nmi1\n%.4f\n',final_nmi11);
% 
% final_acc12=sum(acc_12)/iter;
% final_nmi12=sum(nmi_12)/iter;
% fprintf('fcm_acc2\n%.4f\n',final_acc12);
% fprintf('fcm_nmi2\n%.4f\n',final_nmi12);