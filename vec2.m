vecs = cell(1,length(prods));
e1 = [1 0 0 0];
e2 = [0 1 0 0];
e3 = [0 0 1 0];
e4 = [0 0 0 1];
for k=1:length(prods)
    for i=1:length(prods{k})
        vecs{k}{i} = zeros(1,4);
        for j=1:length(prods{k}{i})
            p = prods{k}{i}{j};
            vecs{k}{i} = vecs{k}{i}+((p==8)-(p==4))*e1+((p==1)-(p==5))*e2+((p==2)-(p==6))*e3+((p==3)-(p==7))*e4;
        end
        
    end
end