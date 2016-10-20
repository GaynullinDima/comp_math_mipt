function [ sol, iter ] = Newton ( func, dfunc, acc, x0, a, b, max_num_iter, multiplicity )
% Newton's method from point x0 in [a,b], a < b.
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

    for i = 1:max_num_iter
        iter = i;
        
        % Check d=0;
        if abs ( dfunc(x) ) < eps
            newx = ( interval(1) + interval(2) ) / 2.0;
            dx = newx - x;
            x = newx;
        else
            dx = -func(x) / dfunc(x) * multiplicity;
            x = x + dx;
        end
        
        % Check pulse.
        if x == history(1)
            x = ( history(1) + history(2) ) / 2.0;
        end
        
        % Update history.
        history(1) = history(2);
        history(2) = x;
        
        % Update interval.
        if ( dfunc(x) > 0 ) && ( func(x) < 0 )
            interval(1) = x;
        elseif ( dfunc(x) > 0 ) && ( func(x) > 0 )
            interval(2) = x;
        elseif ( dfunc(x) < 0 ) && ( func(x) < 0 )
            interval(2) = x;
        elseif (dfunc(x) < 0 ) && ( func(x) > 0 )
            interval(1) = x;
        end
        
        if abs(dx) < acc
            sol = x;
            return;
        end
    end
    
    sol = NaN;
    error ( 'Too many iterations' );
end

