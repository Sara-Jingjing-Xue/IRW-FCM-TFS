function [F, obj1, obj_nfcm, iter,iter_f] = nfcm(U, r, X)
%  X d*n input matrix
%  U n*c initial membership matrix
%  r     fuzzifier parameter

    [~,n] = size(X);
    c=size(U,2);
    F=sparse(U);
    maxiter=30;
    XX=sparse(X'*X);
    X=sparse(X);
    XXd=sparse(diag(XX));
    
    obj_nfcm=zeros(1,maxiter);
    %% initial objective function value     
    G=F.^r;
    center = (G'*X')./(sum(G',2)*ones(1,size(X',2))); % new center
    dist = distfcm(full(center), full(X'));           % fill the distance matrix
    obj_nfcm(1) = sum(sum((dist.^2).*G'));            % objective function
    
for iter=1:maxiter
    tic
    %% s
    s=zeros(c,1);
    G=F.^r;
    sumG=sum(G);
    Z = X*G;
    dd=zeros(1,c);
    for j=1:c
        dd(j)=(Z(:,j)'*Z(:,j))^(1/2);
    end
        
    for j=1:c
        s(j)=dd(j) /sumG(j);
    end       
    %% F
    H=repmat(s',n,1);
    SS=sparse(H.^2);
    F=sparse(F);
    [F,~,obj1,iter_f(iter)] = ff(SS,s,F,XXd,r,100,X);
   %% obj
    G=F.^r;     
    center = G'*X'./(sum(G',2)*ones(1,size(X',2))); % new center
    dist = distfcm(full(center), full(X'));         % fill the distance matrix
    obj_nfcm(iter+1) = sum(sum((dist.^2).*G'));     % objective function
    
   if (iter>1 && abs(obj_nfcm(iter+1)-obj_nfcm(iter)) < 10^-5)
       break;
   end 
   
end
end
    
    
