function X = prod2(prod)
b = cell(1,8);
for k = 1:8 %vpa
    b(k) = {[1+sqrt(2) exp(1i*k*pi/4)*sqrt(2)*sqrt(sqrt(2)+1) ; exp(-1i*k*pi/4)*sqrt(2)*sqrt(sqrt(2)+1) 1+sqrt(2)]};
end
X =b{prod{1}};
for n =2:(length(prod))
X = X * b{prod{n}};
end
end