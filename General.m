b = cell(1,8);
g=2;
for k=1:4*g
    b{k} = [1+sqrt(2) exp(1i*k*pi/2/g)*sqrt(2)*sqrt(sqrt(2)+1); exp(-1i*k*pi/2/g)*sqrt(2)*sqrt(sqrt(2)+1) 11+sqrt(2)]; 
end

% x1 = -2 - sqrt(3);
% y1 = 1 + sqrt(3);
% z1 =  3 + 2*sqrt(3);
% w1 = -1 + sqrt(3);
% f = (3+2*sqrt(3))^(1/4);
% G = cell(1,8);
% G{1} = 1/2*[ -x1 + y1*f z1 + w1*f; -z1 + w1*f   -x1 - y1*f ];
% G{2} = 1/2*[ x1 - w1*f z1 + y1*f; -z1 + y1*f   x1 + w1*f ];
% % for n=2:8
% G{n} = G{1};
% G{n}(1,1) = G{n}(1,1);
% end