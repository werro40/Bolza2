function f = OrbitPlotter2(X)
A = gam(X(1),X(2),X(3),X(4),X(5),X(6));
r = 1/sqrt(2+2*sqrt(2));
c =sqrt((3+2*sqrt(2))/(2+2*sqrt(2)));
theta = linspace(3*pi/4,5*pi/4,100);
octa = cell(1,8);
for j=1:8
    octa(j) = {exp(1i*(j-1)*pi/4)*(c+r*exp(1i*theta))};
end

if(abs(imag(A(1,1)))>1e-3)
    f = {@(x,y) x.^2 + y.^2 - 2/imag(A(1,1)) * ( real(A(1,2)) * y - imag(A(1,2)) * x)+1};
else
    f = {@(x,y) ( real(A(1,2)) * y - imag(A(1,2)) * x)};
end
% oct = cell(1,8);
% for k = 0:7
% oct(k+1) ={@(x,y) x.^2 + y.^2 - 2*x*sqrt(5/4)*cos(k*pi/4) - 2*y*sqrt(5/4)*sin(k*pi/4) + 1};
% end
circ = @(x,y) x.^2 + y.^2 - 1;
fimplicit(circ,'k')
hold on;
% for m = 1:8
%     fimplicit(oct{m}, [-1 1 -1 1],'r');
% end
for m = 1:8
    plot(real(octa{m}),imag(octa{m}),'r')
end
fimplicit(f,[-1 1 -1 1],'b');
end
