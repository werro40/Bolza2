for n1=1:length(vecs)
    disp(n1)
    data = [];
    for j=1:length(vecs{n1})
        for k=1:length(G)
            data = cat(1,data,(vecs{n1}{j}*G{k}));
        end
    end
    data = unique(data,'rows');
    file = sprintf('data/%d.dat',n1);
    dlmwrite(file,data)
end
