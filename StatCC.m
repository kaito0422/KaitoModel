MAX_DISTANCE = 1024;


N0 = sum(l2rdh0);
N1 = sum(l2rdh1);

p1 = N1 / N0;       % 
p2 = N0 / N1;       % 

ccrdh0 = zeros(MAX_DISTANCE, 1);
ccrdh1 = zeros(MAX_DISTANCE, 1);
for rd = 1 : 1024
    rd_tmp0 = round((1 + N1 / N0) * rd);
    rd_tmp1 = round((1 + N0 / N1) * rd);
    if rd_tmp0 > 1024
        ccrdh0(1024) = ccrdh0(1024) + l2rdh0(rd);
    else
        ccrdh0(rd_tmp0) = ccrdh0(rd_tmp0) + l2rdh0(rd);
    end
    
    if rd_tmp1 > 1024
        ccrdh1(1024) = ccrdh1(1024) + l2rdh1(rd);
    else
        ccrdh1(rd_tmp1) = ccrdh1(rd_tmp1) + l2rdh1(rd);
    end
end

plot(ccrdh0);
hold on;
plot(ccrdh1);

CCRDH = zeros(MAX_DISTANCE, 1);
for i = 1 : 1024
    CCRDH(i) = ccrdh0(i) + ccrdh1(i);
end

plot(CCRDH);
hold on;
plot(l2rdh);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






MAX_DISTANCE = 1024;


L2_ASSOC = 16;

% pr(i)表示重用距离大于等于i的概率
Pr = [];
for i = 1:MAX_DISTANCE
    tmp = 0;
    for j = i:MAX_DISTANCE
        tmp = tmp + l2rdh(j);
    end
    Pr(i) = tmp / sum(l2rdh(1:MAX_DISTANCE));
end

% Es(r)表示重用距离为r时，堆栈距离的数学期望
Es = [];
for r = 1:MAX_DISTANCE
    tmp = 0;
    for i = 1:r-1
        tmp = tmp + Pr(i);
    end
    Es(r) = tmp;
end

new_sdh = [];
for s = 1:MAX_DISTANCE
    tmp = 0;
    for r = 1:MAX_DISTANCE
        if round(Es(r)) == s - 1
            tmp = tmp + l2rdh(r);
        end
    end
    new_sdh(s) = tmp;
end

miss_rate_real = sum(new_sdh(L2_ASSOC + 1:MAX_DISTANCE)) / sum(new_sdh);














L2_ASSOC = 16;

% pr(i)表示重用距离大于等于i的概率
Pr = [];
for i = 1:MAX_DISTANCE
    tmp = 0;
    for j = i:MAX_DISTANCE
        tmp = tmp + CCRDH(j);
    end
    Pr(i) = tmp / sum(CCRDH(1:MAX_DISTANCE));
end

% Es(r)表示重用距离为r时，堆栈距离的数学期望
Es = [];
for r = 1:MAX_DISTANCE
    tmp = 0;
    for i = 1:r-1
        tmp = tmp + Pr(i);
    end
    Es(r) = tmp;
end

new_sdh = [];
for s = 1:MAX_DISTANCE
    tmp = 0;
    for r = 1:MAX_DISTANCE
        if round(Es(r)) == s - 1
            tmp = tmp + CCRDH(r);
        end
    end
    new_sdh(s) = tmp;
end

miss_rate_kaito = sum(new_sdh(L2_ASSOC + 1:MAX_DISTANCE)) / sum(new_sdh);








