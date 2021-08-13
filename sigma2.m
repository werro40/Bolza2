function sigma = sigma2(X)
emp = @(n,p) (mod(ceil(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(floor(sqrt(2)*n),2) * floor(sqrt(2)*n)) * mod(p,2)...
    + (mod(floor(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(ceil(sqrt(2)*n),2) * floor(sqrt(2)*n))*mod(p+1,2);

B1 = X(3)*sqrt(2) + emp(X(3),X(5));
B2 = X(4)*sqrt(2) + emp(X(4),X(6));

sigma  = B1*B1 +5*B2*B2 -4*B1*B2;
end