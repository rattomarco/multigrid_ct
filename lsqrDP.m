function x = lsqrDP(A,b,tau,delta,k)
[m,n] = size(A);
U = zeros(m,k+1);
V = zeros(n,k);
B = zeros(k+1,k);

beta = norm(b);
nrmb = beta;
U(:,1) = b/beta;


for j = 1:k
    r =A'*U(:,j);
    if j>1
        r = r - beta*V(:,j-1);
        r = r - V(:,1:j-1)*(V(:,1:j-1)'*r);
    end
    
    alpha = norm(r);
    V(:,j) = r/alpha;
    B(j,j) = alpha;
    p = A*V(:,j);
    p = p-alpha*U(:,j);
    p = p - U(:,1:j)*(U(:,1:j)'*p);
    beta = norm(p);
    B(j+1,j) = beta;
    %%%%%%%%% lsqr %%%%%%%%%
    e = zeros(j+1,1);
    e(1) = nrmb;
    y = B(1:j+1,1:j)\e;
    if norm(B(1:j+1,1:j)*y - e) <= tau*delta || j == k
         x = V(:,1:j)*y;
         break
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    U(:,j+1)=p/beta;
end

