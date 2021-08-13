function Y = conju2(X,k)
Y = zeros(1,6);
emp = @(n,p) (mod(ceil(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(floor(sqrt(2)*n),2) * floor(sqrt(2)*n)) * mod(p,2)...
    + (mod(floor(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(ceil(sqrt(2)*n),2) * floor(sqrt(2)*n))*mod(p+1,2);


Y(1) = X(1);
if(mod(k,2)==1)
    Y(2) = round(5 * X(2) + 4 * emp(X(2),2) + 2*(cos((k-1)*pi/4)+sin((k-1)*pi/4))*( X(3)+ emp( X(3),X(5) ) ) ...
        - 2 * ( cos((k-1)*pi/4) - sin((k-1)*pi/4) ) * ( X(4)+ emp( X(4),X(6) ) ));
    Y(3) = round(3 * X(3) + 2 * emp(X(3),X(5)) -2 * ( (X(3) + emp(X(3),X(5)))*cos(k*pi/2) + sin(k*pi/2)*( X(4) + emp(X(4),X(6)) ) ) + (cos((k-1)*pi/4) + sin((k-1)*pi/4)) * (6 * X(2) +4 * emp(X(2),2)));
    Y(4) = round(3 * X(4) + 2 * emp(X(4),X(6)) -2 * ( (X(3) + emp(X(3),X(5)))*sin(k*pi/2) - cos(k*pi/2)*( X(4) + emp(X(4),X(6)) ) ) - (cos((k-1)*pi/4) - sin((k-1)*pi/4))* (6 * X(2) +4 * emp(X(2),2)));
elseif(mod(k,2)==0)
    Y(2) = round(5 * X(2) + 4 * emp(X(2),2) + 4 * X(3)*sin(k*pi/4)+ 2 * emp(X(3),X(5)) *sin(k*pi/4) ...
        -4*X(4)*cos(k*pi/4) -2*emp(X(4),X(6))*cos(k*pi/4));
    Y(3) = round(3 * X(3) + 2 * emp(X(3),X(5)) -2 * ( ( X(3) + emp(X(3),X(5)) )*cos(k*pi/2) + sin(k*pi/2)*( X(4) + emp(X(4),X(6)) ) ) + sin(k*pi/4)* (8 * X(2) +6 * emp(X(2),2)));
    Y(4) = round(3 * X(4) + 2 * emp(X(4),X(6)) -2 * ( (X(3) + emp(X(3),X(5)))*sin(k*pi/2) - cos(k*pi/2)*( X(4) + emp(X(4),X(6)) ) ) - cos(k*pi/4)* (8 * X(2) +6 * emp(X(2),2)));
end
Y(5) = X(5);
Y(6) = X(6);

end