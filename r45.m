function Y=r45(X,r)
em = @(n) mod(ceil(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(floor(sqrt(2)*n),2) * floor(sqrt(2)*n); %must be odd
em2 = @(n) mod(floor(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(ceil(sqrt(2)*n),2) * floor(sqrt(2)*n); %must be even
for i=1:r
        Xp5 = mod(X(3) - X(4),2);
        Xp6 = mod(X(3) + X(4),2);
        Xp3 = ( em(X(3)) * X(5) + em2(X(3)) * (1-X(5)) )/2 - (em(X(4)) * X(5) + em2(X(4)) * (1-X(5)))/2;
        Xp4 = ( em(X(3)) * X(6) + em2(X(3)) * (1-X(6)) )/2 + (em(X(4)) * X(6) + em2(X(4)) * (1-X(6)))/2;
        X(3) = Xp3;
        X(4) = Xp4;
        X(5) = Xp5;
        X(6) = Xp6;
end
Y=X;
end