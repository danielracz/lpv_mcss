function PlotData(filename)

load(filename);

figure()
plot(NumberOfTs(3:end),ge_errs_1(3:end),NumberOfTs(3:end), Bounds_1(3:end)+emp_errs_1(3:end))
legend('True error','PAC bound+empirical error')
xlabel('N')
title('L1 loss')

figure()
plot(NumberOfTs,ge_errs_1-emp_errs_1,NumberOfTs, Bounds_1)
legend('Generalization gap','PAC bound')
xlabel('N')
title('L1 loss')

figure()
plot(NumberOfTs,ge_errs-emp_errs,NumberOfTs, Bounds_2)
legend('Generalization gap','PAC bound')
title('L2 loss"')
xlabel('N')


figure()
plot(NumberOfTs,ge_errs,NumberOfTs, Bounds_2+emp_errs)
legend('True error','PAC bound+empirical error')
title('L2 loss')
xlabel('N')


end