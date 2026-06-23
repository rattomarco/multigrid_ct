function [P,R] = projectors(xsize,bsize,lev,method)

P = cell(lev+1,1);
R = cell(lev+1,1);


if method == 1  % [0 1 1]
    for i = 0:lev
        % R
        b1 = ceil(bsize(1)/(2^i));
        b2 = ceil(bsize(2)/(2^i));
        Cr=spdiags([zeros(b1,1) ones(b1,1) ones(b1,1)]/2, -1:1,b1,b1);
        Cc=spdiags([zeros(b2,1) ones(b2,1) ones(b2,1)]/2, -1:1,b2,b2);
        Kr=speye(b1);
        Kr=Kr(1:2:end,:);
        Kc=speye(b2);
        Kc=Kc(1:2:end,:);
        K=kron(Kc,Kr);
        R{i+1}=K*(kron(Cc,Cr)/(2^i));

        % P
        x1 = ceil(xsize(1)/(2^i));
        x2 = ceil(xsize(2)/(2^i));
        Cr=spdiags([zeros(x1,1) ones(x1,1) ones(x1,1)]/2, -1:1,x1,x1);
        Cc=spdiags([zeros(x2,1) ones(x2,1) ones(x2,1)]/2, -1:1,x2,x2);
        K=speye(x1);
        K=K(1:2:end,:);
        K=kron(K,K);
        P{i+1}=(K*(kron(Cc,Cr)/(2^i)))';
    end
end

if method == 2
    for i = 0:lev
        % R
        b1 = ceil(bsize(1)/(2^i));
        b2 = ceil(bsize(2)/(2^i));
        Cr=spdiags([ones(b1,1) 2*ones(b1,1) ones(b1,1)]/4, -1:1,b1,b1);
        Cc=spdiags([ones(b2,1) 2*ones(b2,1) ones(b2,1)]/4, -1:1,b2,b2);
        Kr=speye(b1);
        Kr=Kr(1:2:end,:);
        Kc=speye(b2);
        Kc=Kc(1:2:end,:);
        K=kron(Kc,Kr);
        R{i+1}=K*(kron(Cc,Cr)/(2^i));

        % P
        x1 = ceil(xsize(1)/(2^i));
        x2 = ceil(xsize(2)/(2^i));
        Cr=spdiags([ones(x1,1) 2*ones(x1,1) ones(x1,1)]/4, -1:1,x1,x1);
        Cc=spdiags([ones(x2,1) 2*ones(x2,1) ones(x2,1)]/4, -1:1,x2,x2);
        K=speye(x1);
        K=K(1:2:end,:);
        K=kron(K,K);
        P{i+1}=(K*(kron(Cc,Cr)/(2^i)))';
    end
end

if method == 3 % [0 1 3 3 1]
    for i = 0:lev
        % R
        b1 = ceil(bsize(1)/(2^i));
        b2 = ceil(bsize(2)/(2^i));
        Cr=spdiags([zeros(b1,1) ones(b1,1) 3*ones(b1,1) 3*ones(b1,1) ones(b1,1)]/8, -2:2,b1,b1);
        Cc=spdiags([zeros(b2,1) ones(b2,1) 3*ones(b2,1) 3*ones(b2,1) ones(b2,1)]/8, -2:2,b2,b2);
        Kr=speye(b1);
        Kr=Kr(1:2:end,:);
        Kc=speye(b2);
        Kc=Kc(1:2:end,:);
        K=kron(Kc,Kr);
        R{i+1}=K*(kron(Cc,Cr)/(2^i));

        % P
        x1 = ceil(xsize(1)/(2^i));
        x2 = ceil(xsize(2)/(2^i));
        Cr=spdiags([zeros(x1,1) ones(x1,1) 3*ones(x1,1) 3*ones(x1,1) ones(x1,1)]/8, -2:2,x1,x1);
        Cc=spdiags([zeros(x1,1) ones(x1,1) 3*ones(x1,1) 3*ones(x1,1) ones(x1,1)]/8, -2:2,x2,x2);
        K=speye(x1);
        K=K(1:2:end,:);
        K=kron(K,K);
        P{i+1}=(K*(kron(Cc,Cr)/(2^i)))';
    end
end

if method == 4 % [1 4 6 4 1]
    for i = 0:lev
        % R
        b1 = ceil(bsize(1)/(2^i));
        b2 = ceil(bsize(2)/(2^i));
        Cr=spdiags([ones(b1,1) 4*ones(b1,1) 6*ones(b1,1) 4*ones(b1,1) ones(b1,1)]/16, -2:2,b1,b1);
        Cc=spdiags([ones(b2,1) 4*ones(b2,1) 6*ones(b2,1) 4*ones(b2,1) ones(b2,1)]/16, -2:2,b2,b2);
        Kr=speye(b1);
        Kr=Kr(1:2:end,:);
        Kc=speye(b2);
        Kc=Kc(1:2:end,:);
        K=kron(Kc,Kr);
        R{i+1}=K*(kron(Cc,Cr)/(2^i));

        % P
        x1 = ceil(xsize(1)/(2^i));
        x2 = ceil(xsize(2)/(2^i));
        Cr=spdiags([ones(x1,1) 4*ones(x1,1) 6*ones(x1,1) 4*ones(x1,1) ones(x1,1)]/16, -2:2,x1,x1);
        Cc=spdiags([ones(x1,1) 4*ones(x1,1) 6*ones(x1,1) 4*ones(x1,1) ones(x1,1)]/16, -2:2,x2,x2);
        K=speye(x1);
        K=K(1:2:end,:);
        K=kron(K,K);
        P{i+1}=(K*(kron(Cc,Cr)/(2^i)))';
    end
end

if method == 10
    P1 = diag(ones(xsize(1),1)) + diag(ones(xsize(1)-1,1),-1);
    P1 = P1(:,1:2:end);
    P2 = kron(P1,P1);
    P2=sparse(P2);

    R=speye(size(A,1));
    R=R(1:4:end,:);
end