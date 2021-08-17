function X = Identifier(M)
A1 = real(M(1,1));
A2 = imag(M(1,1));
B1 = real(M(1,2))/sqrt(sqrt(2)-1);
B2 = imag(M(1,2))/sqrt(sqrt(2)-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load EvenLookup1.mat EvenLookup1
load OddLookup1.mat OddLookup1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%        figure out n1:
n1p1 = find(abs(A1 - OddLookup1)<1e-3);
if(isempty(n1p1))
    n1p2 = find(abs(A1 + OddLookup1)<1e-3);
    if(~isempty(n1p2))
        n1 = 1-(n1p2);
    else
        n1 = 42.7;
    end
else
    n1 = n1p1-1;
end

%          figure out n2:
n2p1 = find(abs(A2 - EvenLookup1)<1e-3);

if(isempty(n2p1))
    n2p2 = find(abs(A2 + EvenLookup1)<1e-3);
    if(~isempty(n2p2))
        n2 = 1-(n2p2);
    else
        n2 = 42.7;
    end
else
    n2 = (n2p1)-1;
end

%          figure out n3:
n3p1 = find(abs(B1 - OddLookup1)<1e-2);
n3p2 = find(abs(B1 + OddLookup1)<1e-2);
n3p3 = find(abs(B1 - EvenLookup1)<1e-2);
n3p4 = find(abs(B1 + EvenLookup1)<1e-2);
if(~isempty(n3p1))
    n3 = (n3p1)-1;
    parity1 = 1;
elseif(~isempty(n3p2))
    n3 = 1-(n3p2);
    parity1 =1;
elseif(~isempty(n3p3))
    n3 = (n3p3)-1;
    parity1 =0;
elseif(~isempty(n3p4))
    n3 = 1-(n3p4);
    parity1 =0;
else
    n3 =42.7;
    parity1 = 42.7;
end

%          figure out n4:
n4p1 = find(abs(B2 - OddLookup1)<1e-3);
n4p2 = find(abs(B2 + OddLookup1)<1e-3);
n4p3 = find(abs(B2 - EvenLookup1)<1e-3);
n4p4 = find(abs(B2 + EvenLookup1)<1e-3);
if(~isempty(n4p1))
    n4 = (n4p1)-1;
    parity2 = 1;
elseif(~isempty(n4p2))
    n4 = 1-(n4p2);
    parity2 =1;
elseif(~isempty(n4p3))
    n4 = (n4p3)-1;
    parity2 =0;
elseif(~isempty(n4p4))
    n4 = 1-(n4p4);
    parity2 =0;
else
    n4 =42.7;
    parity2 =42.7;
end

if(parity1~=parity2)
    if((n3~=0))
        parity2 = parity1;
    elseif((n4~=0))
        parity1 = parity2;
    elseif((n3 ==0)&& (n4==0))
        parity1 =0;
        parity2 =0;
    end
end


X = [n1 n2 n3 n4 parity1 parity2];
end