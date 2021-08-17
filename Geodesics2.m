% function g = Geodesics2(N) %#codegen

pick = 1;
N=80;
prods = cell(1,N);
Length = zeros(N,1);
P =[1 0 0 0 0 0; 0 -1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0;0 0 0 0 0 1];
% Q = 5e2;
g = zeros(N,1); %multiplicities
s5000= zeros(N,1);
errProd = [];
% trouble = 0;
cl=0;
% em = @(n) mod(ceil(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(floor(sqrt(2)*n),2) * floor(sqrt(2)*n); %must be odd
% em2 = @(n) mod(floor(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(ceil(sqrt(2)*n),2) * floor(sqrt(2)*n); %must be even
emp = @(n,p) (mod(ceil(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(floor(sqrt(2)*n),2) * floor(sqrt(2)*n)) * mod(p,2)...
    + (mod(floor(sqrt(2)*n),2) * ceil(sqrt(2)*n) + mod(ceil(sqrt(2)*n),2) * floor(sqrt(2)*n))*mod(p+1,2);
d = @(A,B) (A-B)*(A-B)';


t=sqrt(sqrt(2)-1);
% load('LengthBolza.mat');
BadGeos2 = cell(1,16);
BadGeos2{13} = [1 -4 6 3 0 0];
BadGeos2{14} = [ 1     4     6     2     1     1];
BadGeos2{15} = [ 1     4     6     3     0     0];
BadGeos2{16} = [ 1    -4     6     2     1     1];
for n1 = 1:N
    XXX = cell(0);

    XX = cell(0);
    disp(n1)
    
    RR = cell(1);
    A1 = n1*sqrt(2) + emp(n1,1) ; %new
    
    Length(n1) = 2 * acosh( A1 );%new %vpa
    special = logical(abs(round(Length(n1)/Length(1))- Length(n1)/Length(1))<1e-8);
    
    
    
    delta = ( 1 + (sqrt(2)) )*( 1 + (sqrt(2)) )*( 1 + (sqrt(2)) ) * ( A1*A1 - 1 );%new %vpa
    indexI = 0;
    for index =0:1 %were missing all the B2 at p=1 term
        
        n3 = 0;
        
        while (n3 <(3*sqrt(delta)-1)/sqrt(8))
            
            
            %disp("B1:")
            B1 = n3*sqrt(2) + emp(n3,index) ;%new
            n4 = floor(1/4*max(0, real(sqrt(A1*A1/t/t-B1*B1-1/t/t)-1)));
            while n4<=n3
                
                %disp("B2:")
                B2 =  n4*sqrt(2) + emp(n4,index) ;%new
                
                sigma  = B1 * B1 + 5 * B2 * B2 - 4 * B1 * B2;
                
                if logical( sigma <= delta-1*(-1).^(special)*1e-12) %is it in1side the octagon
                    for neg =1:2
                        A2 = (-1).^neg*sqrt(1 + (sqrt(2)-1) * (B1*B1 +B2*B2) - A1*A1);
                        
                        
                        %                         Valid = sum(abs(A2 - EvenLookup1)<1e-6); %nice, must agree within a millionth
                        %                         ind = find(abs(A2 - EvenLookup1)<1e-6);
                        %                         A2 = EvenLookup1(ind);
                        %                         n2 = ind-1;
                        guess = floor((A2+1)/sqrt(8));
                        n2 = guess;
                        if(isreal(A2))
                            if((abs(A2-guess*sqrt(2)-emp(guess,0))<1e-5))
                                
                                Valid1 = 1;
                            else
                                Valid1= 0;
                            end
                        else
                            Valid1= 0;
                        end
                        X =[n1 n2 n3 n4 index index];
                        %                         M = [ A1 + 1i * A2 sqrt(sqrt(2)-1) * ( B1 + 1i * B2); sqrt(sqrt(2)-1) * (B1 - 1i * B2) A1 - 1i * A2];
                        if(Valid1&&(d(BadGeos2{13},X)>1e-8)&&(d(BadGeos2{14},X)>1e-8)&&(d(BadGeos2{15},X)>1e-8)&&(d(BadGeos2{16},X)>1e-8)||((~special)&&Valid1)) %is it a valid geodesic
                            
                            %which # of valid matrix is this
                            fin = 0; % is the series ended
                            
%                             M = cell(1);
                            X = cell(1);
                            
                            B = cell(1,8); %conjugated
                            R = cell(1,8); %conjugate and rotate
%                             MR = cell(1);
                            XR = cell(1);
                            
                            Conditions = zeros(1,8);
                            
                            X{1} =[n1 n2 n3 n4 index index];
                            XR{1} = X{1};
                            if(d(X{1},[ 1 -1 2 1 1 1])<1e-2)
                                 disp("o")
                            end
                            if(~isempty(XX))
                                Mrot =cell(1,8);
                                for  r= 1:8
                                    Mrot{r} =r45(X{1},r);
                                end
                                for m=1:(length(XX))
                                    for r = 1:8
                                        if(d(XX{m},Mrot{r})<1e-8)
                                            fin = 1;
                                            cl = 1;
                                        end
                                    end
                                    if(d(r45(XX{m}*P,4),X{1})<1e-8)
                                        fin = 1;
                                        cl=1;
                                    end

                                end
                            end
                            if(fin==1)
                                if(cl==1)
                                    prods{n1}(end)=[];
                                    XXX(end) =[];
                                    cl=0;
                                end
                                break
                            end
                            
                            indexI = indexI + 1; %which # of valid matrix is this
                           
                            %                             MM(indexI) = M(1);
                            XX{indexI} = X{1};
                            % Making sure this was not used before
                            
                            for k=1:8
                                B{k} = conju2(X{1},k);
                                
                            end
                            i=0; %do any elements match M{1}
                            for k=1:8
                                R{k} = rot2(B{k});
                                
                                
                                if(sigma2(R{k})-delta<=-1*(-1).^(special)*1e-9)
                                    if((d2(BadGeos2{13},R{k})>1e-8)&&(d2(BadGeos2{14},R{k})>1e-8)&&(d2(BadGeos2{15},R{k})>1e-8)&&(d2(BadGeos2{16},R{k})>1e-8)||(~special))
                                        i=i+1;
                                        Conditions(i) = k;
                                        
                                    end
                                end
                                
                            end
                            % % % %            5%%%%%%%%%%%%%%%%%%%%%%%%
                            r=1;
                            if(Conditions(1)==0)
                                
                                for k=1:8
                                    
                                    if(sigma2(R{k})-delta<=1)
                                        
                                        i=i+1;
                                        Conditions(i) = k;
                                        
                                    end
                                    
                                end
                                %                                 Conditions(1)=1;
                                %                                 R(1) = X(1);
                                %prods{n1}{end+1}{1} = Conditions(1);
                            end
%                             possFin =0;
%                             for f=1:8
%                                 if(d(B{f},X{1}))
%                                     possFin =1;
%                                 end
%                             end
                            if(Conditions(2)==0)
                                pick = 1;
                                X{r+1} = B{Conditions(pick)};
                                XR{r} = R{Conditions(pick)};
                                prods{n1}{end+1}{1} = Conditions(pick);
                                XXX{end+1}{1} = X{pick};
                            elseif((sigma2(R{Conditions(2)})-delta==0)&&special)
                                pick = 2;
                                X{r+1} = B{Conditions(pick)};
                                XR{r} = R{Conditions(pick)};
                                prods{n1}{end+1}{1} = Conditions(pick);
                                XXX{end+1}{1} = X{pick};
                            else
                                pick = 1;
                                X{r+1} = B{Conditions(pick)};
                                XR{r} = R{Conditions(pick)};
                                prods{n1}{end+1}{1} = Conditions(pick);
                                XXX{end+1}{1} = X{pick};
                            end
                            
                            
                            %check if M{2} == M{1}
                            if(d(R{Conditions(1)},X{1})< 1e-8)
                                fin = 1;
%                                 disp("Original:")
%                                 disp(X{1})
%                                 disp("Predicted:")
%                                 disp(Identifier(prod2(prods{end})))
%                                 if(trouble)
%                                     
%                                     disp(X{1});
%                                     figure
%                                     for i=1:length(XR)
%                                         OrbitPlotter2(X{i});
%                                     end
%                                 end
                                errProd = [errProd logical( d2( abs(X{1}), abs(Identifier(prod2(prods{n1}{end}))) ) <1e-2)];
                                if(errProd(end)==0)
                                    disp("Error")
                                end
                                RR{indexI} = X{1};
                                XR{1}= X{1};
                                
                                g(n1) = g(n1) + 8;
                                
                            else
                                X{2} = B{Conditions(1)};
                                XR{1} = R{Conditions(1)};
                                RR{indexI} = R{Conditions(1)};
                                %prods{end}{end+1} = Conditions(1);
                                 XXX{end}{2} = X{2};
                            end
                            if(length(XX)>1)
                                for z =1:(length(XX)-1)
                                    if((d(XX{z},X{1})<1e-8)||(d(RR{z},XR{1})<1e-8))
                                        fin = 1;
                                        cl=1;
                                    end
                                end
                            end
                            % % % %            5%%%%%%%%%%%%%%%%%%%%%%%%
                            if(fin==1)
                                if(cl==1)
                                    prods{n1}(end)=[];
                                    XXX(end) =[];
                                    cl=0;
                                end
                                break
%                                 prods(end)=[];
%                                 XXX(end) =[];
                            end
                            r=2;
                            
                            while(fin == 0)
                                i=0;
                                Conditions = zeros(1,8);
                                
                                
                                
                                for k=1:8
                                    B{k} = conju2(X{r},k);
                                    
                                end
                                
                                for k=1:8
                                    R{k} = rot2(B{k});
                                    
                                    
                                    if(sigma2(R{k})-delta<=-1e-9*(-1).^(special))
                                        if((d2(BadGeos2{13},R{k})>1e-9)&&(d2(BadGeos2{14},R{k})>1e-9)&&(d2(BadGeos2{15},R{k})>1e-9)&&(d2(BadGeos2{16},R{k})>1e-9)||(~special))
                                            if(d(B{k},X{r-1})>1e-3)
                                                i=i+1;
                                                Conditions(i) = k;
                                            end
                                            
                                        end
                                    end
                                    
                                end
                                
                                if((Conditions(1)==0))
                                    
                                    for k=1:8
                                        
                                        if(sigma2(R{k})-delta<=1)
                                            
                                            
                                                i=i+1;
                                                Conditions(i) = k;
                                            
                                            
                                        end
                                        
                                    end
                                    
                                    %                                     X(r+1) = M(1);
                                    %                                     XR(r) = X(r);
                                    %prods{end}{end+1} = 0;
                                    %                                     g(n1) = g(n1) +4;
                                    %                                     break
                                end
                                possFin =zeros(1,8);
                                for f=1:8
                                    if(d(B{f},X{1})<1e-3)
                                        possFin(f) =1;
                                        
                                    end
                                end
                                possBad =zeros(1,8);
                                for f=1:8
                                    if(d(B{f},X{r-1})<1e-3)
                                        possBad(f) =1;
                                        
                                    end
                                end
                                Cond = Conditions(Conditions>0);
                                score =zeros(1,length(Cond));
                                for h=1:length(Cond)
                                    if((possFin(Cond(h))==0)&&(possBad(Cond(h))==0))
                                        score(h) = 3*(~special)+4*(special);
                                    elseif((possFin(Cond(h))==1)&&(possBad(Cond(h))==0))
                                        score(h) = 4*(~special)+3*(special);
                                    elseif((possFin(Cond(h))==1)&&(possBad(Cond(h))==1))
                                        score(h) = 2;
                                    elseif((possFin(Cond(h))==0)&&(possBad(Cond(h))==1))
                                        score(h) = 1;
                                    end
                                end
                                maxes = find(score==max(score));
                                pick = maxes(1);
                                sortScore = sort(score);
                                X{r+1} = B{Conditions(pick)};
                                XR{r} = R{Conditions(pick)};
                                prods{n1}{end}{end+1} = Conditions(pick);
                                XXX{end}{end+1} = X{r+1};
%                                 if((length(Cond)>1)&&(abs(prods{end}{end}-prods{end}{end-1})==4))
%                                     prods{end}{end} = Conditions(sortScore(2));
%                                 end


                                if(length(XXX)>1)
                                    for z =1:(length(XX)-1)
                                        if((d(r45(XX{z}*P,4),X{r})<1e-8)||(d(XX{z},X{r})<1e-8)||(d(RR{z},XR{r})<1e-8))
                                            fin =1;
                                            cl=1;
                                        end
                                    end
                                    Mrot =cell(1,8);
                                    for  s1= 7:9
                                        Mrot(s1-6) ={r45(X{r},s1)};
                                    end
                                    
                                    for z =1:(length(XXX)-1)
                                        for s1=1:3
                                            for y =1:length(XXX{z})
                                                if((d(r45(XXX{z}{y}*P,4),Mrot{s1})<1e-8)||(d(XXX{z}{y},Mrot{s1})<1e-8))
                                                    fin =1;
                                                    cl=1;
                                                end
                                            end
                                        end
                                    end
                                    if(fin==1)
                                        prods{n1}(end) = [];
                                        XXX(end) = [];
                                        break
                                    end
                                end

                                
                                
                                %%%%%%%%%%%%
                                if(d(B{Conditions(pick)},X{1})< 1e-3)
                                     fin = 1;
                                     
                                     errProd = [errProd logical( d2( abs(X{1}), abs(Identifier(prod2(prods{n1}{end}))) ) <1e-2)];
                                     if(errProd(end)==0)
                                         disp("Error")
                                     end
                                     %                                     disp("Original:")
%                                     disp(X{1})
%                                     disp("Predicted:")
%                                     disp(Identifier(prod2(prods{end})))
%                                     if(trouble)
%                                         disp(X{1});
%                                         figure
%                                         for i=1:length(XR)
%                                             OrbitPlotter2(X{i});
%                                             
%                                         end
%                                     end
                                    Nb2 = 0;
                                    Nb4 = 0;
                                    Nb8 = 0;
                                    if(mod(length(XR),2) == 0)
                                        %Nb2 = logical( (d(X{1},XR{length(XR)/2+1}) < 1e-3)||(d(r45(X{1},1),XR{length(XR)/2+1}) < 1e-3)||(d(r45(X{1},7),XR{length(XR)/2+1}) < 1e-3))*8;
                                        Nb2 = logical( (d(X{1},XR{length(XR)/2}) < 1e-3)||(d(X{1},r45(XR{length(XR)/2},1)) < 1e-3)||(d(X{1},r45(XR{length(XR)/2},7)) < 1e-3))*4;
                                    end
                                    if((mod(length(XR),4) == 0))
                                        %Nb4 = logical( (d(X{1},XR{length(XR)/4+1}) < 1e-3)||(d(r45(X{1},1),XR{length(XR)/4+1}) < 1e-3)||(d(r45(X{1},7),XR{length(XR)/4+1}) < 1e-3))*4;
                                        Nb2 = logical( (d(X{1},XR{length(XR)/4}) < 1e-3))*2;
                                    end
                                    if((mod(length(XR),8) == 0))
                                        %Nb8 = logical( (d(X{1},XR{length(XR)/8+1}) < 1e-3)||(d(r45(X{1},1),XR{length(XR)/8+1}) < 1e-3)||(d(r45(X{1},7),XR{length(XR)/8+1}) < 1e-3))*2;
                                        Nb2 = logical( (d(X{1},XR{length(XR)/8}) < 1e-3))*1;
                                    end
                                    if(Nb8~=0)
                                        Nb = 1;
                                    elseif(Nb4~=0)
                                        Nb = 2;
                                    elseif(Nb2~=0)
                                        Nb = 4;
                                    else
                                        Nb = 8;
                                    end
                                    
                                    g(n1) = g(n1) + Nb;
                                    
                                else
                                    RR{indexI} = R{Conditions(pick)};
                                end
                                
                                
                                r = r+1;
                                if(r>35)
                                    fin = 1;
                                    
                                    Nb = 4;
                                    g(n1) =g(n1) + Nb;
                                    s5000(n1) = s5000(n1)+1;
                                end
                                
                            end
                            
                            
                        end
                    end
                end
                n4 = n4 + 1; %new
            end
            n3 = n3 + 1; %new
        end
    end
    
end
for n=1:N
    for k=2:5
        for j =1:n
            if(abs(k*Length(j)-Length(n))<1e-3)
                g(n) = g(n) - g(j);
            end
        end
    end
end
%g([4,24,48,72])=0;
zeroList =[4,24,48,72,140,160,176,184,200,224,288,432,456,472,496,704,728,816,856,880,1024,1088,1112,1128,1136,1152,1360,1384,1400,1424,1488];
zero = zeroList(zeroList<=N);
g(zero)=0;
%disp("Deviation:")
%disp(d(g',gacc(1:N)'))
%disp(find(abs(g-gacc(1:N))>0)')
% end
%g([4,24,48,72,140,160,176,184,200,224,288,432,456,472,496,704,728,816,856,880,1024,1088,1112,1128,1136,1152,1360,1384,1400,1424,1488]) = 0;