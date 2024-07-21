plot(log2(data(:,1)), log2(data(:,2)), 'o');
hold on;
plot(log2(data(:,1)), log2(data(:,3)), 's');
hold on;
plot(log2(data(:,1)), log2(data(:,4)), '^');
hold on;
plot(log2(data(:,1)), log2(data(:,5)), 'v');
hold on;
plot(log2(data(:,1)), log2(data(:,6)), '>');
hold on;
plot(log2(data(:,1)), log2(data(:,7)), '<');
hold on;

legend("field3_1dp", "field3_1dp_t", "field3_3dp", "field3_r", "field3_map", "field3_mdspan", 'Interpreter', 'none');
xlabel("log2(scale)");
ylabel("log2(time)(s)");
title("performance test of different field3 implementations");