function [x0, Ymean, Ahat, Bhat, Chat]=dytex(Y,n,nv)
tau = size(Y,2);
Ymean = mean(Y,2);

[U,S,V] = svd(Y-Ymean*ones(1,tau),0);

Chat = U(:,1:n);
Xhat = S(1:n,1:n)*V(:,1:n)';
x0=Xhat(:,1);
Ahat = Xhat(:,2:tau)*pinv(Xhat(:,1:(tau-1)));
% Ahat(Ahat>1)=1;
% Ahat(Ahat<-1)=-1;
Vhat = Xhat(:,2:tau)-Ahat*Xhat(:,1:(tau-1));
[Uv,Sv,Vv] = svd(Vhat,0);
Bhat = Uv(:,1:nv)*Sv(1:nv,1:nv)/sqrt(tau-1);
end