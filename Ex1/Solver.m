function sol = Solver ( A, b )
    [ Alu, ordr ] = LUdec_modified ( A );
    
    sol = LUsol_modified ( Alu, b, ordr );
    
end


