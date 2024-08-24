close all ;clear;
num = 100; %行走次数
time = 100; %实验次数
total1 = zeros(time,1);
total2 = zeros(time,1);
total3 = zeros(time,1);
total4 = zeros(time,1);
final1 = zeros(time,1);
final2 = zeros(time,1);
final3 = zeros(time,1);
final4 = zeros(time,1);
for i = 1:time
    [RMSE1,RMSE2,RMSE3,RMSE4,variance1,variance2,variance3,variance4] = trail1();
    total1(i) = RMSE1;
    total2(i) = RMSE2;
    total3(i) = RMSE3;
    total4(i) = RMSE4;
    final1(i) = variance1;
    final2(i) = variance2;
    final3(i) = variance3;
    final4(i) = variance4;
end
figure;
plot(total1,'r-*','linewidth',1.4);
xlabel('实验次数'); ylabel('均方根误差/m');
grid on;
hold on;
plot(total2,'g-*');
hold on;
plot(total3,'b-*');
hold on;
plot(total4,'y-*');
hold on;
save('RMSE','total1','total2','total3','total4','final1','final2','final3','final4');
for i = 1:time
%     plot([i,i],[total1(i)-final1(i)/10,total1(i)+final1(i)/10],'m-d');
%     hold on;
    plot([i,i],[total2(i)-final2(i)/10,total2(i)+final2(i)],'b-d');
    hold on;
end
legend('未重采样算法', '误差椭圆重采样','分层重采样', '系统重采样');
hold off;
figure;
cdfplot(total1);
hold on;
cdfplot(total2);
xlabel('Distance Error'); ylabel('CDF');
hold off;