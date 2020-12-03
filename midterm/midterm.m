figure(1)
syn = struct('t1',1,'gsyn',1*1e-4,'loc',0.06,'vsyn',-20);
[~, v1_1] = mid1(syn);
syn = struct('t1',3,'gsyn',1*1e-4,'loc',0.04,'vsyn',70);
[t, v1_2] = mid1(syn);
v1(1,:) = v1_1;
v1(2,:) = v1_2;
prt(t,v1)
figure(2)
syn = struct('t1',3,'gsyn',1*1e-4,'loc',0.04,'vsyn',-20);
[~, v1_1] = mid1(syn);
syn = struct('t1',1,'gsyn',1*1e-4,'loc',0.06,'vsyn',70);
[t, v1_2] = mid1(syn);
v1(1,:) = v1_1;
v1(2,:) = v1_2;
prt(t,v1)