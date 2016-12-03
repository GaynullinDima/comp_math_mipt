function x = LUsol_modified( A, b, ordr )
% Solver with LU-method.
    if size ( b, 2 ) > 1
        % Convert to column.
        b = b';
    end
    
    b = ChangeStringsByOrder ( b, ordr );
    
    n = length ( b );
    
    for k = 2:n
        b(k) = b(k) - A(k,1:k-1) * b(1:k-1);
    end
    
    for k = n:-1:1
        b(k) = ( b(k) - A(k,k+1:n) * b(k+1:n) ) / A(k,k);
    end
    
    x = b;
end

function B = ChangeStringsByOrder ( A, ordr )
% Changes strings to theirs order.
    n = size ( A, 1 );
    m = size ( A, 2 );
    B = zeros ( n, m );
    
    for i = 1:n
        cur = ordr(i);
        B(i,:) = A(cur,:);
    end
end
