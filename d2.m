function x = d2(A,B)
emp = @(n,p) (mod(ceil(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(floor(sqrt(2)*n),2) * floor(sqrt(2)*n)) * mod(p,2)...
    + (mod(floor(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(ceil(sqrt(2)*n),2) * floor(sqrt(2)*n))*mod(p+1,2);
if(A(2)==0||B(2)==0)
    x= abs(A(2)-B(2)) + abs(A(3)-B(3)) + abs(A(4)-B(4)) + abs(A(5)-B(5)) + abs(A(6)-B(6));
else
    A2a = A(2)*sqrt(2) + emp(A(2),2);
    A2b = B(2)*sqrt(2) + emp(B(2),2); 
    B1a = A(3)*sqrt(2) + emp(A(3),A(5));
    B1b = B(3)*sqrt(2) + emp(B(3),B(5));
    B2a = A(4)*sqrt(2) + emp(A(4),A(6));
    B2b = B(4)*sqrt(2) + emp(B(4),B(6));
    
    x = (B1a/A2a - B1b/A2b)^2 + (B2a/A2a - B2b/A2b)^2 ; 
end



end