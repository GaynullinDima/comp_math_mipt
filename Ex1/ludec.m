function A = ludec( A)
%[L/U] form
n = size(A, 1);
for k = 1:n-1
    for i = k+1:n
        if A(i,k)~=0
            lambda = A(i,k)/A(k,k);
            A(i,k+1:n) = A(i,k+1:n) - lambda*A(k,k+1:n);
            A(i,k) = lambda;
        end
    end
end
end


