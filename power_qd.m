close all
qd8=VarName3(85:3960,:);
qd3=VarName1(85:3960,:);
qd2=VarName2(85:3960,:);


qd8=qd8/qd8(3500);
qd3=qd3/qd3(3500);
qd2=qd2/qd2(3500);

plot(qd8);
hold on
plot(qd3);
plot(qd2);

set(gca, 'YScale', 'log')
axis([0 3876 0.5 100])