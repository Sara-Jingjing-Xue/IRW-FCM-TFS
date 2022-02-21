
function [F,A,obj1,i] = ff(SS,s,F,XXd,r,maxiter,X )

    [~,c]=size(F);
    obj1=zeros(1,maxiter);
    dd=zeros(1,c);
    G=F.^r;
       
    for i=1:maxiter             
        Z = X*G;
        B=X'*Z;
        for j=1:c
            dd(j)=(Z(:,j)'*Z(:,j))^(1/2);
        end    
        D=diag(1./dd);
        A=B*D;  % a_j
        
        dist = repmat(XXd,1,c) + SS - 2*A*diag(s);
        tmp = dist.^(-1/(r-1)); 
        F= tmp./(sum(tmp,2)*ones(1, c));        
%%
        G = F.^r;    
        sumG=sum(G,2); 
        obj1(i) =sumG'*XXd+ sum(G)*s.^2 - 2*s'*diag(G'*A); 

        if (i>1 && abs(obj1(i)-obj1(i-1)) < 10^-5)
           break;
        end 
    end
end