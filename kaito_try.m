MAX_DISTANCE = 1024;
L1_ASSOC = 2;

prs = zeros(MAX_DISTANCE, MAX_DISTANCE);
for r = 1 : MAX_DISTANCE
    sum_tmp = 0;
    for s = 1 : MAX_DISTANCE
        sum_tmp = sum_tmp + rst(r, s);
    end
    for s = 1 : MAX_DISTANCE
        if sum_tmp == 0
            prs(r, s) = 0;
        else
            prs(r, s) = rst(r, s) / sum_tmp;
        end
    end
end

miss_rdh = zeros(MAX_DISTANCE, 1);
for r = 1 : MAX_DISTANCE
    prop = 0;
    for i = 1 : L1_ASSOC
        prop = prop + prs(r, i);
    end
    miss_rdh(r) = drdh(r) * (1 - prop);
end

% Phit_num(i, j)表示在L1重用距离为i的reference中，对应的reuse epoch中有多少个references在
% L1 cache中是hit的
Phit_num = zeros(MAX_DISTANCE, MAX_DISTANCE);
for rd = 1 : MAX_DISTANCE
    sum_tmp = 0;
    for n = 1 : rd
        sum_tmp = sum_tmp + hit_rdh(rd, n);
    end
    for n = 1 : rd
        if sum_tmp == 0
            Phit_num(rd, n) = 0;
        else
            Phit_num(rd, n) = hit_rdh(rd, n) / sum_tmp;
        end
    end
end

kaito_l2rdh = zeros(MAX_DISTANCE, 1);
for rd = 1 : MAX_DISTANCE
    for i = 1 : rd              % 最后在L2的重用距离为i，则需要在L1有rd - i个命中
        kaito_l2rdh(i) = kaito_l2rdh(i) + miss_rdh(rd) * Phit_num(rd, rd - i + 1);
    end
end

plot(l2rdh);
hold on;
plot(kaito_l2rdh);


plot(l2rdh);
hold on;
plot(miss_rdh);

%kaito_l2rdh(MAX_DISTANCE) = kaito_l2rdh(MAX_DISTANCE) + miss_rdh(MAX_DISTANCE) * Phit_num(MAX_DISTANCE, 1);

% plot(l2rdh);
% hold on;
% plot(kaito_l2rdh);

sum1 = sum(l2rdh);
sum2 = sum(miss_rdh);
sum3 = sum(kaito_l2rdh);

% plot(miss_rdh);
% hold on;
% plot(kaito_l2rdh);



plot(l2rdh);
hold on;
plot(kaito_l2rdh);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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

% new_l2rdh = zeros(MAX_DISTANCE, 1);
% for rd = 1 : MAX_DISTANCE - 1
%     if mod(rd, 2) == 0
%         new_l2rdh(rd / 2) = new_l2rdh(rd / 2) + kaito_l2rdh(rd);
%     else
%         new_l2rdh((rd + 1) / 2) = new_l2rdh((rd + 1) / 2) + kaito_l2rdh(rd);
%     end
% end

Psame = 0.5; % L1 set num / L2 set num

new_l2rdh = zeros(MAX_DISTANCE, 1);
for rd1 = 0 : MAX_DISTANCE - 2
    for rd2 = 0 : rd1
        new_l2rdh(rd2 + 1) = new_l2rdh(rd2 + 1) + Cmn(rd1 + 1, rd2 + 1) * ((Psame)^rd2) * ((1 - Psame)^(rd1 - rd2)) * kaito_l2rdh(rd1 + 1);
    end
end
%new_l2rdh(MAX_DISTANCE) = kaito_l2rdh(MAX_DISTANCE) / 4;

plot(l2rdh);
hold on;
plot(new_l2rdh);

% plot(l2rdh);
% hold on;
% plot(kaito_l2rdh);

sum1 = sum(l2rdh);
sum2 = sum(new_l2rdh);

plot(drdh);
plot(miss_rdh);



















