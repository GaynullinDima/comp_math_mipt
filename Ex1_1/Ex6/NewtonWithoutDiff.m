function [ sol, iter ] = NewtonWithoutDiff ( func, acc, x0, a, b, max_num_iter, multiplicity )
% Newton's method without given differential from point x0 in [a,b], a < b.
% Done:
%       1) Pulse.
%       2) d=0;
%       3) Multiplicity of roots.
    
    eps = acc/1000.0;

    x = x0;
    
    % Interval [a, b] - for d=0.
    interval = 1:2;
    interval(1) = a;
    interval(2) = b;
    
    % History - to protect from pulse = биение.
    history = 1:2;
    history(1) = x0;
    history(2) = x0;
    
    % Diff
    u = 1:6;
    h = acc;
    l = 2;
    m = 3;

    for iter = 1:max_num_iter
        % Count d.
        
        for n = -l:m
            u(n+l+1) = func(x+h*n);
        end
        [ d, mu ] = get_diff ( u, l, m, h );
        
        % Update interval.
        if ( d > 0 ) && ( func(x) < 0 )
            interval(1) = x;
        elseif ( d > 0 ) && ( func(x) > 0 )
            interval(2) = x;
        elseif ( d < 0 ) && ( func(x) < 0 )
            interval(2) = x;
        elseif ( d < 0 ) && ( func(x) > 0 )
            interval(1) = x;
        end
        
        % Check d=0;
        
        if abs ( d ) < eps
            newx = ( interval(1) + interval(2) ) / 2.0;
            dx = newx - x;
            x = newx;
        else
            dx = -func(x) / d * multiplicity;
            x = x + dx;
        end
        
        % Check pulse.
        if x == history(1)
            x = ( history(1) + history(2) ) / 2.0;
        end
        
        % Update history.
        history(1) = history(2);
        history(2) = x;
        
        
        
        if abs(dx) < acc
            sol = x;
            return;
        end
    end
    
    sol = NaN;
    error ( 'Too many iterations' );
end

function [ diff, mu ] = get_diff ( u, l, m, h )
    n = length ( u );
    v = -l:m;
    
    A = fliplr(vander(v))';
    
    b = zeros ( n, 1 );
    b(2) = 1;
    
    alpha = A\b;
    
    diff = 1 / h * sum(dot(alpha,u));
    
    mu = cond ( A );
end
