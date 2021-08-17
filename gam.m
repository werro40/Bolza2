function M = gam(n1,n2,n3,n4,e1,e2)
    em = @(n) mod(ceil(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(floor(sqrt(2)*n),2) * floor(sqrt(2)*n); %must be odd
    em2 = @(n) mod(floor(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(ceil(sqrt(2)*n),2) * floor(sqrt(2)*n);
    A1 = em(n1) +sqrt(2)*n1;
    A2 = em2(n2) +sqrt(2)*n2;
    if(e1==0)
        B1 = em2(n3) +sqrt(2)*n3;
        
    else
        B1 = em(n3) +sqrt(2)*n3;
        
    end
    
    if(e2==0)
        
        B2 = em2(n4) +sqrt(2)*n4;
    else
        
        B2 = em(n4) +sqrt(2)*n4;
    end

    
    M = [ A1 + 1i * A2 sqrt(sqrt(2)-1) * ( B1 + 1i * B2); sqrt(sqrt(2)-1) * (B1 - 1i * B2) A1 - 1i * A2];
end