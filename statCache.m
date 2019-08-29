N = sum(l2rdh);
L2_ASSOC = 16;

LEFT = [];
RIGHT = [];
i = 1;

error = 10000;
l2_miss_rate = 0;
for R = 0.00001:0.00001:0.99999
    left = N * R;
    right = 0;
    for r = 1:1024
        f = 1 - (1 - 1 / L2_ASSOC)^(r*R);
        right = right + l2rdh(r)*f;
    end
    LEFT(i) = left;
    RIGHT(i) = right;
    i = i + 1;
    
    if abs(left - right) < error
        error = abs(left - right);
        l2_miss_rate = R;
    end
end





N = sum(new_l2rdh);
L2_ASSOC = 16;

LEFT = [];
RIGHT = [];
i = 1;

error = 10000;
l2_miss_rate = 0;
for R = 0.00001:0.00001:0.99999
    left = N * R;
    right = 0;
    for r = 1:1024
        f = 1 - (1 - 1 / L2_ASSOC)^(r*R);
        right = right + new_l2rdh(r)*f;
    end
    LEFT(i) = left;
    RIGHT(i) = right;
    i = i + 1;
    
    if abs(left - right) < error
        error = abs(left - right);
        l2_miss_rate = R;
    end
end

