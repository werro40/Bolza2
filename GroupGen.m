R = [0 0 0 -1;1 0 0 0;0 1 0 0; 0 0 1 0];
S = [0 1 0 0;1 0 0 0;0 0 0 -1; 0 0 -1 0];
T = [0 -1 1 -1;-1 0 1 -1;-1 1 0 -1; -1 1 -1 0];
U = [0 -1 0 0;1 -1 0 0;1 0 -1 1; 0 1 -1 0];
G = cell(0);
dM = @(A,B) trace((A-B)*(A-B)');
for r=0:7
    for s=0:1
        for t=0:1
            for u=0:2
                M = R^r*S^s*T^t*U^u;
                rep =0;
                for i=1:length(G)
                    if(dM(M,G{i})<1e-3)
                        rep = rep + 1;
                    end
                end
                if(rep==0)
                    G{end+1} = M;
                end
            end
        end
    end
end