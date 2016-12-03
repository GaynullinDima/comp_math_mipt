function [ A, ordr ] = LUdec_modified ( A )
    % [L/U] form with choosing max element in column.
    % ordr - order of strings.
    
    n = size ( A, 1 );
    
    ordr = 1:n;
    
    for k = 1:n-1
        % Find string "str" with main element 
        % - max in column lower than string #k.
        str = MaxLowElement ( A, k );
        
        if str ~= k
            % Change strings "str" and "k".
            A = ChangeStrings ( A, k, str );

            % Change ordr for returning original order of strings later.
            buffer = ordr(k);
            ordr(k) = ordr(str);
            ordr(str) = buffer;
        end
        
        % LU of this string.
        
        for i = k+1:n
            if A(i,k) ~= 0
                lambda = A(i,k) / A(k,k);
                
                A(i,k+1:n) = A(i,k+1:n) - lambda * A(k,k+1:n);
                A(i,k) = lambda;
            end
        end
    end
end

function num_str = MaxLowElement ( A, cur_str )
% Finds Max element in column cur_str after string cur_str-1.
    sz = size ( A, 1 );

    cur_max = A(cur_str,cur_str);
    num_str = cur_str;
    
    for i = cur_str+1:sz
        if A(i,cur_str) > cur_max
            cur_max = A(i,cur_str);
            num_str = i;
        end
    end
end

function B = ChangeStrings ( A, cur_str, max_str )
% Changes strings cur_str and max_str in A.
    n = size ( A, 1 );

    B = A;
    B(cur_str,1:n) = A(max_str,1:n);
    B(max_str,1:n) = A(cur_str,1:n);
end


