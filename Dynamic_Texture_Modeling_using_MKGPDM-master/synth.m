function I=synth(x0,Ymean,Ahat, Bhat, Chat,tau)
[n,nv]=size(Bhat);
X(:,1) = x0;
I=zeros(size(Ymean,1),tau);
for t=1:tau
    X(:,t+1)=Ahat * X(:,t)+Bhat*randn(nv,1);
    I(:,t)=Chat*X(:,t)+Ymean;
    %I(:,t) = (I(:,t) - floor(min(I(:))))./(ceil(max(I(:)))-floor(min(I(:))));
end
end