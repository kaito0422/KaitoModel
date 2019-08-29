

MAX_DISTANCE = 1024;
L2_ASSOC = 8;

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

miss_rate = sum(new_sdh(L2_ASSOC + 1:MAX_DISTANCE)) / sum(new_sdh);












% 
% 
% 
% MAX_DISTANCE = 1024;
% L2_ASSOC = 16;
% 
% % pr(i)表示重用距离大于等于i的概率
% Pr = [];
% for i = 1:MAX_DISTANCE
%     tmp = 0;
%     for j = i:MAX_DISTANCE
%         tmp = tmp + new_l2rdh(j);
%     end
%     Pr(i) = tmp / sum(new_l2rdh(1:1024));
% end
% 
% % Es(r)表示重用距离为r时，堆栈距离的数学期望
% Es = [];
% for r = 1:MAX_DISTANCE
%     tmp = 0;
%     for i = 1:r-1
%         tmp = tmp + Pr(i);
%     end
%     Es(r) = tmp;
% end
% 
% new_sdh = [];
% for s = 1:MAX_DISTANCE
%     tmp = 0;
%     for r = 1:MAX_DISTANCE
%         if round(Es(r)) == s - 1
%             tmp = tmp + new_l2rdh(r);
%         end
%     end
%     new_sdh(s) = tmp;
% end
% 
% miss_rate = sum(new_sdh(L2_ASSOC + 1:MAX_DISTANCE)) / sum(new_sdh);
% 

















MAX_DISTANCE = 1024;
L2_ASSOC = 8;

% pr(i)表示重用距离大于等于i的概率
Pr = [];
for i = 1:MAX_DISTANCE
    tmp = 0;
    for j = i:MAX_DISTANCE
        tmp = tmp + kaito_l2rdh(j);
    end
    Pr(i) = tmp / sum(kaito_l2rdh(1:MAX_DISTANCE));
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
            tmp = tmp + kaito_l2rdh(r);
        end
    end
    new_sdh(s) = tmp;
end

miss_rate_lru = sum(new_sdh(L2_ASSOC + 1:MAX_DISTANCE)) / sum(new_sdh);








