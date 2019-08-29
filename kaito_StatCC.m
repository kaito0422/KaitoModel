MAX_DISTANCE = 1024;
L1_ASSOC = 4;
L2_ASSOC = 16;

prs0 = zeros(MAX_DISTANCE, MAX_DISTANCE);
for r = 1 : MAX_DISTANCE
    sum_tmp = 0;
    for s = 1 : MAX_DISTANCE
        sum_tmp = sum_tmp + rst0(r, s);
    end
    for s = 1 : MAX_DISTANCE
        if sum_tmp == 0
            prs0(r, s) = 0;
        else
            prs0(r, s) = rst0(r, s) / sum_tmp;
        end
    end
end

miss_rdh0 = zeros(MAX_DISTANCE, 1);
for r = 1 : MAX_DISTANCE
    prop = 0;
    for i = 1 : L1_ASSOC
        prop = prop + prs0(r, i);
    end
    miss_rdh0(r) = drdh0(r) * (1 - prop);
end

% Phit_num(i, j)表示在L1重用距离为i的reference中，对应的reuse epoch中有多少个references在
% L1 cache中是hit的
Phit_num0 = zeros(MAX_DISTANCE, MAX_DISTANCE);
for rd = 1 : MAX_DISTANCE
    sum_tmp = 0;
    for n = 1 : rd
        sum_tmp = sum_tmp + hit_rdh0(rd, n);
    end
    for n = 1 : rd
        if sum_tmp == 0
            Phit_num0(rd, n) = 0;
        else
            Phit_num0(rd, n) = hit_rdh0(rd, n) / sum_tmp;
        end
    end
end

kaito_l2rdh0 = zeros(MAX_DISTANCE, 1);
for rd = 1 : MAX_DISTANCE
    for i = 1 : rd              % 最后在L2的重用距离为i，则需要在L1有rd - i个命中
        kaito_l2rdh0(i) = kaito_l2rdh0(i) + miss_rdh0(rd) * Phit_num0(rd, rd - i + 1);
    end
end

%kaito_l2rdh(MAX_DISTANCE) = kaito_l2rdh(MAX_DISTANCE) + miss_rdh(MAX_DISTANCE) * Phit_num(MAX_DISTANCE, 1);

plot(kaito_l2rdh0);
hold on;
plot(l2rdh0);

Cmn = zeros(MAX_DISTANCE, MAX_DISTANCE);
for m = 1 : MAX_DISTANCE
    for n = 1 : m
        if n == 1 || n == m
            Cmn(m, n) = 1;
        else
            Cmn(m, n) = Cmn(m - 1, n - 1) + Cmn(m - 1, n);
        end
    end
end

Psame = 0.5; % L1 set num / L2 set num

new_l2rdh0 = zeros(MAX_DISTANCE, 1);
for rd1 = 0 : MAX_DISTANCE - 2
    for rd2 = 0 : rd1
        new_l2rdh0(rd2 + 1) = new_l2rdh0(rd2 + 1) + Cmn(rd1 + 1, rd2 + 1) * ((Psame)^rd2) * ((1 - Psame)^(rd1 - rd2)) * kaito_l2rdh0(rd1 + 1);
    end
end

plot(new_l2rdh0);
hold on;
plot(l2rdh0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


prs1 = zeros(MAX_DISTANCE, MAX_DISTANCE);
for r = 1 : MAX_DISTANCE
    sum_tmp = 0;
    for s = 1 : MAX_DISTANCE
        sum_tmp = sum_tmp + rst1(r, s);
    end
    for s = 1 : MAX_DISTANCE
        if sum_tmp == 0
            prs1(r, s) = 0;
        else
            prs1(r, s) = rst1(r, s) / sum_tmp;
        end
    end
end

miss_rdh1 = zeros(MAX_DISTANCE, 1);
for r = 1 : MAX_DISTANCE
    prop = 0;
    for i = 1 : L1_ASSOC
        prop = prop + prs1(r, i);
    end
    miss_rdh1(r) = drdh0(r) * (1 - prop);
end

% Phit_num(i, j)表示在L1重用距离为i的reference中，对应的reuse epoch中有多少个references在
% L1 cache中是hit的
Phit_num1 = zeros(MAX_DISTANCE, MAX_DISTANCE);
for rd = 1 : MAX_DISTANCE
    sum_tmp = 0;
    for n = 1 : rd
        sum_tmp = sum_tmp + hit_rdh1(rd, n);
    end
    for n = 1 : rd
        if sum_tmp == 0
            Phit_num1(rd, n) = 0;
        else
            Phit_num1(rd, n) = hit_rdh1(rd, n) / sum_tmp;
        end
    end
end

kaito_l2rdh1 = zeros(MAX_DISTANCE, 1);
for rd = 1 : MAX_DISTANCE
    for i = 1 : rd              % 最后在L2的重用距离为i，则需要在L1有rd - i个命中
        kaito_l2rdh1(i) = kaito_l2rdh1(i) + miss_rdh1(rd) * Phit_num1(rd, rd - i + 1);
    end
end

%kaito_l2rdh(MAX_DISTANCE) = kaito_l2rdh(MAX_DISTANCE) + miss_rdh(MAX_DISTANCE) * Phit_num(MAX_DISTANCE, 1);

plot(kaito_l2rdh1);
hold on;
plot(l2rdh1);

Cmn = zeros(MAX_DISTANCE, MAX_DISTANCE);
for m = 1 : MAX_DISTANCE
    for n = 1 : m
        if n == 1 || n == m
            Cmn(m, n) = 1;
        else
            Cmn(m, n) = Cmn(m - 1, n - 1) + Cmn(m - 1, n);
        end
    end
end

Psame = 0.5; % L1 set num / L2 set num

new_l2rdh1 = zeros(MAX_DISTANCE, 1);
for rd1 = 0 : MAX_DISTANCE - 2
    for rd2 = 0 : rd1
        new_l2rdh1(rd2 + 1) = new_l2rdh1(rd2 + 1) + Cmn(rd1 + 1, rd2 + 1) * ((Psame)^rd2) * ((1 - Psame)^(rd1 - rd2)) * kaito_l2rdh1(rd1 + 1);
    end
end

plot(new_l2rdh1);
hold on;
plot(l2rdh1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(new_l2rdh0);
hold on;
plot(new_l2rdh1);

plot(l2rdh0);
hold on;
plot(l2rdh1);

N0 = sum(new_l2rdh0);
N1 = sum(new_l2rdh1);

N0 = sum(kaito_l2rdh0);
N1 = sum(kaito_l2rdh1);

p1 = N1 / N0;       % 
p2 = N0 / N1;       % 

ccrdh0 = zeros(MAX_DISTANCE, 1);
ccrdh1 = zeros(MAX_DISTANCE, 1);
for rd = 1 : 1024
    rd_tmp0 = round((1 + N1 / N0) * rd);
    rd_tmp1 = round((1 + N0 / N1) * rd);
    if rd_tmp0 > 1024
        ccrdh0(1024) = ccrdh0(1024) + kaito_l2rdh0(rd);
    else
        ccrdh0(rd_tmp0) = ccrdh0(rd_tmp0) + kaito_l2rdh0(rd);
    end
    
    if rd_tmp1 > 1024
        ccrdh1(1024) = ccrdh1(1024) + kaito_l2rdh1(rd);
    else
        ccrdh1(rd_tmp1) = ccrdh1(rd_tmp1) + kaito_l2rdh1(rd);
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

sum1 = sum(CCRDH);
sum2 = sum(l2rdh);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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







