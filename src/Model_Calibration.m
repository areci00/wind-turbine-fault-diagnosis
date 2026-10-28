clear 
close all
clc

load  WT_Fault_free_data.mat


t=Omega_r_m1(:,1);

%% plot data

figure (1)
subplot (5,3,1)
hold on
plot (t,Omega_r_m1(:,2))
title ('Omega_rm1')
subplot (5,3,2)
hold on
plot (t,Omega_r_m2(:,2))
title ('Omega_rm2')
subplot (5,3,3)
hold on
plot (t,Omega_g_m1(:,2))
title ('Omega_gm1')
subplot (5,3,4)
hold on
plot (t,Omega_g_m2(:,2))
title ('Omega_gm2')
subplot (5,3,5)
hold on
plot (t,Beta_1_m1(:,2))
title ('Beta_1m1')
subplot (5,3,6)
hold on
plot (t,Beta_1_m2(:,2))
title ('Beta_1m2')
subplot (5,3,7)
hold on
plot (t,Beta_2_m1(:,2))
title ('Beta_2m1')
subplot (5,3,8)
hold on
plot (t,Beta_2_m2(:,2))
title ('Beta_2m2')
subplot (5,3,9)
hold on
plot (t,Beta_3_m1(:,2))
title ('Beta_3m1')
subplot (5,3,10)
hold on
plot (t,Beta_3_m2(:,2))
title ('Beta_3m2')
subplot (5,3,11)
hold on
plot (t,Tau_g_m(:,2))
title ('Tau_gm')
subplot (5,3,12)
hold on
plot (t,P_g_m(:,2))
title ('P_gm')

% Control signals: T_g_r, Tau_r, Beta_r
subplot (5,3,13)
hold on
plot (t,T_g_r(:,2))
title ('T_g_r')

subplot (5,3,14)
hold on
plot (t,T_g_r(:,2))
title ('Tau_r')

subplot (5,3,15)
hold on
plot (t,Beta_r(:,2))
title ('Beta_r')


%% Filter

Filter=1; %0=>No filter 
fc=0.1;
if Filter==1
    
    n_butt=2;
    Ts=t(2)-t(1);
    [n_butt,d_butt]=butter(2,fc);
    filter_data=tf(n_butt,d_butt,Ts);
    Beta_1_m1(:,2)=lsim (filter_data,Beta_1_m1(:,2));
    Beta_1_m2(:,2)=lsim (filter_data,Beta_1_m2(:,2));
    Beta_2_m1(:,2)=lsim (filter_data,Beta_2_m1(:,2));
    Beta_2_m2(:,2)=lsim (filter_data,Beta_2_m2(:,2));
    Beta_3_m1(:,2)=lsim (filter_data,Beta_3_m1(:,2));
    Beta_3_m2(:,2)=lsim (filter_data,Beta_3_m2(:,2));
    Omega_r_m1(:,2)=lsim (filter_data,Omega_r_m1(:,2));
    Omega_r_m2(:,2)=lsim (filter_data,Omega_r_m2(:,2));
    Omega_g_m1(:,2)=lsim (filter_data,Omega_g_m1(:,2));
    Omega_g_m2(:,2)=lsim (filter_data,Omega_g_m2(:,2));
    Tau_g_m(:,2)=lsim (filter_data,Tau_g_m(:,2));
    P_g_m(:,2)=lsim (filter_data,P_g_m(:,2));
    T_g_r(:,2)=lsim (filter_data,T_g_r(:,2));
    Beta_r(:,2)=lsim (filter_data,Beta_r(:,2));
    Tau_r(:,2)=lsim (filter_data,Tau_r(:,2));
    
end
%%

figure (1)
subplot (5,3,1)
plot (t,Omega_r_m1(:,2),'r')
subplot (5,3,2)
plot (t,Omega_r_m2(:,2),'r')
subplot (5,3,3)
plot (t,Omega_g_m1(:,2),'r')
subplot (5,3,4)
plot (t,Omega_g_m2(:,2),'r')
subplot (5,3,5)
plot (t,Beta_1_m1(:,2),'r')
subplot (5,3,6)
plot (t,Beta_1_m2(:,2),'r')
subplot (5,3,7)
plot (t,Beta_2_m1(:,2),'r')
subplot (5,3,8)
plot (t,Beta_2_m2(:,2),'r')
subplot (5,3,9)
plot (t,Beta_3_m1(:,2),'r')
subplot (5,3,10)
plot (t,Beta_3_m2(:,2),'r')
subplot (5,3,11)
plot (t,Tau_g_m(:,2),'r')
subplot (5,3,12)
plot (t,P_g_m(:,2),'r')
subplot (5,3,13)
plot (t,T_g_r(:,2),'r')
subplot (5,3,14)
plot (t,T_g_r(:,2),'r')
subplot (5,3,15)
plot (t,Beta_r(:,2),'r')


%% Generating residuals
beta_th=1.1; %10% extra to the thresholds


%Inicialization variables
N=max(size(Beta_1_m1));
r1=zeros(N,1);
r2=zeros(N,1);
Est_Omega_r=zeros(N,1);
r3=zeros(N,1);
r4=zeros(N,1);
Est_Omega_g=zeros(N,1);
r5=zeros(N,1);
r6=zeros(N,1);
Est_Beta_1=zeros(N,1);
r7=zeros(N,1);
r8=zeros(N,1);
Est_Beta_2=zeros(N,1);
r9=zeros(N,1);
r10=zeros(N,1);
Est_Beta_3=zeros(N,1);
r11=zeros(N,1);
Est_Tau_g=zeros(N,1);
r12=zeros(N,1);

load WT_Parameter_models.mat

T_ini=6000; %transient of filters

% Residual 1 (static)
r1(T_ini:end)=Omega_r_m1(T_ini:end,2)-Omega_r_m2(T_ini:end,2);
figure (101)
subplot (2,1,1)
plot(t,Omega_r_m1(:,2),'b',t,Omega_r_m2(:,2),'r')
subplot (2,1,2)
plot (t,r1)

th_1=beta_th*max(abs(r1(1000:end)))

% Residual 2 (dynamic)
Est_Omega_r(1)=mean(Omega_r_m2(1:100,2)); %Initial conditions

L_2=0; %Observer gain, L_2=(0,a2_1) (simulation/prediction)



for i=2:N,
    Est_Omega_r(i)=a2_1*Est_Omega_r(i-1)+b2_1*Omega_g_m2(i-1,2)+c2_1*Tau_r(i-1,2)+d2_1*Tau_g_m(i-1,2)+ L_2*(Omega_r_m2(i-1,2)-Est_Omega_r(i-1));
end

r2(T_ini:end)=Omega_r_m2(T_ini:end,2)-Est_Omega_r(T_ini:end);

figure (102)
subplot (2,1,1)
plot (t,Omega_r_m2(:,2),'b',t,Est_Omega_r,'r')

subplot (2,1,2)
plot (t,r2)

th_2=beta_th*max(abs(r2))

%% Residual 3 (static)

r3(T_ini:end)=Omega_g_m1(T_ini:end,2)-Omega_g_m2(T_ini:end,2);

figure (103)
subplot (2,1,1)
plot(t,Omega_g_m1(:,2),'b',t,Omega_g_m2(:,2),'r')
subplot (2,1,2)
plot (t,r3)

th_3=beta_th*max(abs(r3))


%% Residual 4  (dynamic)

Est_Omega_g(1)=Omega_g_m2(1,2);

L_4=0*a4_1;  %Observer gain, L_2=(0,a4_1) (simulation/prediction)

for i=2:N,
    Est_Omega_g(i)= a4_1*Est_Omega_g(i-1)+b4_1*Omega_r_m2(i-1,2)+c4_1*Tau_r(i-1,2)+d4_1*Tau_g_m(i-1,2)+L_4*(Omega_g_m2(i-1,2)-Est_Omega_g(i-1));
    
end

r4(T_ini:end)=Omega_g_m2(T_ini:end,2)-Est_Omega_g(T_ini:end);


figure (104)
subplot (2,1,1)
plot (t,Omega_g_m2(:,2),'b',t,Est_Omega_g,'r')

subplot (2,1,2)
plot (t,r4)

th_4=beta_th*max(abs(r4))


%% Residu 5

r5=Beta_1_m1(:,2)-Beta_1_m2(:,2);

figure (105)
subplot (2,1,1)
plot(t,Beta_1_m1(:,2),'b',t,Beta_1_m2(:,2),'r')
subplot (2,1,2)
plot (Beta_1_m1(:,1),r5)

th_5=beta_th*max(abs(r5(6000:end)))


%% Residuo 6

L_6=0*a6_1;  %Observer gain, L_6=(0,a4_1) (simulation/prediction)

for i=3:N,
    Est_Beta_1(i)=a6_1*Est_Beta_1(i-1)+a6_2*Est_Beta_1(i-2)+b6_1*Beta_r(i-1,2)+b6_2*Beta_r(i-2,2)+L_6*(Beta_1_m2(i-1,2)-Est_Beta_1(i-1));
    
end

r6(T_ini:end)=Beta_1_m2(T_ini:end,2)-Est_Beta_1(T_ini:end);

figure (106)
subplot (2,1,1)
plot (t,Beta_1_m2(:,2),'b',t,Est_Beta_1,'r')

subplot (2,1,2)
plot (t,r6)

th_6=beta_th*max(abs(r6(6000:end)))

%% Residual 7

r7=Beta_2_m1(:,2)-Beta_2_m2(:,2);

figure (107)
subplot (2,1,1)
plot(t,Beta_2_m1(:,2),'b',t,Beta_2_m2(:,2),'r')
subplot (2,1,2)
plot (t,r7)

th_7=beta_th*max(abs(r7(6000:end)))


%% Residuo 8
L_8=0*a8_1;
Est_Beta_2(1:2,1)=Beta_2_m2(1:2,2);

for i=3:N,
    Est_Beta_2(i)=a8_1*Est_Beta_2(i-1)+a8_2*Est_Beta_2(i-2)+b8_1*Beta_r(i-1,2)+b8_2*Beta_r(i-2,2)+L_8*(Beta_2_m2(i-1,2)-Est_Beta_2(i-1));
end

r8(T_ini:end)=Beta_2_m2(T_ini:end,2)-Est_Beta_2(T_ini:end);



figure (108)
subplot (2,1,1)
plot (t,Beta_2_m2(:,2),'b',t,Est_Beta_2,'r')

subplot (2,1,2)
plot (t,r8)

th_8=beta_th*max(abs(r8))

%% Residu 9

r9=Beta_3_m1(:,2)-Beta_3_m2(:,2);

figure (109)
subplot (2,1,1)
plot(t,Beta_3_m1(:,2),'b',t,Beta_3_m2(:,2),'r')
subplot (2,1,2)
plot (t,r9)

th_9=beta_th*max(abs(r9(6000:end)))


%% Residual 10

L_10=0*a10_1;
Est_Beta_3(1:2,1)=Beta_3_m2(1:2,2);

for i=3:N,
    Est_Beta_3(i)=a10_1*Est_Beta_3(i-1)+a10_2*Est_Beta_3(i-2)+b10_1*Beta_r(i-1,2)+b10_2*Beta_r(i-2,2)+L_10*(Beta_3_m2(i-1,2)-Est_Beta_3(i-1));
end

r10(T_ini:end)=Beta_3_m2(T_ini:end,2)-Est_Beta_3(T_ini:end);

figure (110)
subplot (2,1,1)
plot (t,Beta_3_m2(:,2),'b',t,Est_Beta_3,'r')

subplot (2,1,2)
plot (t,r10)

th_10=beta_th*max(abs(r10))

%% Residual 11

L_11=0*a11_1;

Est_Tau_g(1)=Tau_g_m(1,2);

for i=2:N,
    Est_Tau_g(i)=a11_1*Est_Tau_g(i-1)+b11_1*T_g_r(i-1,2)+L_11*(Tau_g_m(i-1,2)-Est_Tau_g(i-1));
end

r11(T_ini:end)=Tau_g_m(T_ini:end,2)-Est_Tau_g(T_ini:end);
figure (111)
subplot (2,1,1)
plot (t,Tau_g_m(:,2),'b',t,Est_Tau_g,'r')

subplot (2,1,2)
plot (t,r11)

th_11=beta_th*max(abs(r11))

%% Residual 12

eta=0.98;
Est_P_g(:,1)=eta*Omega_g_m1(:,2).*Tau_g_m(:,2);

r12(T_ini:end)=P_g_m(T_ini:end,2)-Est_P_g(T_ini:end);

figure (112)
subplot (2,1,1)
plot(t,P_g_m(:,2),'b',t,Est_P_g(:,1),'r')
subplot (2,1,2)
plot (t,r12)

th_12=beta_th*max(abs(r12))


figure (200)
subplot (4,3,1)
hold on
plot (t,r1)
plot (t,th_1*ones(size(t)),t,-th_1*ones(size(t)),'r')
title ('r_1')
subplot (4,3,2)
hold on
plot (t,r2)
plot (t,th_2*ones(size(t)),t,-th_2*ones(size(t)),'r')
title ('r_2')
subplot (4,3,3)
hold on
plot (t,r3)
plot (t,th_3*ones(size(t)),t,-th_3*ones(size(t)),'r')
title ('r_3')
subplot (4,3,4)
hold on
plot (t,r4)
plot (t,th_4*ones(size(t)),t,-th_4*ones(size(t)),'r')
title ('r_4')
subplot (4,3,5)
hold on
plot (t,r5)
plot (t,th_5*ones(size(t)),t,-th_5*ones(size(t)),'r')
title ('r_5')
subplot (4,3,6)
hold on
plot (t,r6)
plot (t,th_6*ones(size(t)),t,-th_6*ones(size(t)),'r')
title ('r_6')
subplot (4,3,7)
hold on
plot (t,r7)
plot (t,th_7*ones(size(t)),t,-th_7*ones(size(t)),'r')
title ('r_7')
subplot (4,3,8)
hold on
plot (t,r8)
title ('r_8')
plot (t,th_8*ones(size(t)),t,-th_8*ones(size(t)),'r')
subplot (4,3,9)
hold on
plot (t,r9)
title ('r_9')
plot (t,th_9*ones(size(t)),t,-th_9*ones(size(t)),'r')
subplot (4,3,10)
hold on
plot (t,r10)
plot (t,th_10*ones(size(t)),t,-th_10*ones(size(t)),'r')
title ('r_{10}')
subplot (4,3,11)
hold on
plot (t,r11)
plot (t,th_11*ones(size(t)),t,-th_11*ones(size(t)),'r')
title ('r_{11}')
subplot (4,3,12)
hold on
plot (t,r12)
plot (t,th_12*ones(size(t)),t,-th_12*ones(size(t)),'r')
title ('r_{12}')

save thresholds Filter fc th_1 th_2 th_3 th_4 th_5 th_6 th_7 th_8 th_9 th_10 th_11 th_12 L_2 L_4 L_6 L_8 L_10 L_11 T_ini


 

