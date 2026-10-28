clear
close all
clc



load  WT_Faulty_data.mat


load thresholds % Thresholds, Filter parameters and Observer Gains
t=Omega_r_m1(:,1);

%% plot data

figure (1)
subplot (12,1,1)
plot (t,Omega_r_m1(:,2))
title ('Omega_rm1')
subplot (12,1,2)
plot (t,Omega_r_m2(:,2))
title ('Omega_rm2')
subplot (12,1,3)
plot (t,Omega_g_m1(:,2))
title ('Omega_gm1')
subplot (12,1,4)
plot (t,Omega_g_m2(:,2))
title ('Omega_gm2')
subplot (12,1,5)
plot (t,Beta_1_m1(:,2))
title ('Beta_1m1')
subplot (12,1,6)
plot (t,Beta_1_m2(:,2))
title ('Beta_1m2')
subplot (12,1,7)
plot (t,Beta_2_m1(:,2))
title ('Beta_2m1')
subplot (12,1,8)
plot (t,Beta_2_m2(:,2))
title ('Beta_2m2')
subplot (12,1,9)
plot (t,Beta_3_m1(:,2))
title ('Beta_3m1')
subplot (12,1,10)
plot (t,Beta_3_m2(:,2))
title ('Beta_3m2')
subplot (12,1,11)
plot (t,Tau_g_m(:,2))
title ('Tau_gm')
subplot (12,1,12)
plot (t,P_g_m(:,2))
title ('P_gm')

%% Filter


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
end
%%

%% Generating residuals

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

fi1=r1;
fi2=r2;
fi3=r3;
fi4=r4;
fi5=r5;
fi6=r6;
fi7=r7;
fi8=r8;
fi9=r9;
fi10=r10;
fi11=r11;
fi12=r12;
Fault=fi1;



load WT_Parameter_models.mat

% Residual 1 (static)
r1(T_ini:end)=Omega_r_m1(T_ini:end,2)-Omega_r_m2(T_ini:end,2);
% Residual 2 (dynamic)
Est_Omega_r(1)=mean(Omega_r_m2(1:100,2)); %Initial conditions
for i=2:N,
    Est_Omega_r(i)=a2_1*Est_Omega_r(i-1)+b2_1*Omega_g_m2(i-1,2)+c2_1*Tau_r(i-1,2)+d2_1*Tau_g_m(i-1,2)+ L_2*(Omega_r_m2(i-1,2)-Est_Omega_r(i-1));
end
r2(T_ini:end)=Omega_r_m2(T_ini:end,2)-Est_Omega_r(T_ini:end);
%% Residual 3 (static)
r3(T_ini:end)=Omega_g_m1(T_ini:end,2)-Omega_g_m2(T_ini:end,2);
%% Residual 4  (dynamic)
Est_Omega_g(1)=Omega_g_m2(1,2);
for i=2:N,
    Est_Omega_g(i)= a4_1*Est_Omega_g(i-1)+b4_1*Omega_r_m2(i-1,2)+c4_1*Tau_r(i-1,2)+d4_1*Tau_g_m(i-1,2)+L_4*(Omega_g_m2(i-1,2)-Est_Omega_g(i-1));
end
r4(T_ini:end)=Omega_g_m2(T_ini:end,2)-Est_Omega_g(T_ini:end);
%% Residu 5
r5(T_ini:end)=Beta_1_m1(T_ini:end,2)-Beta_1_m2(T_ini:end,2);

%% Residuo 6
for i=3:N,
    Est_Beta_1(i)=a6_1*Est_Beta_1(i-1)+a6_2*Est_Beta_1(i-2)+b6_1*Beta_r(i-1,2)+b6_2*Beta_r(i-2,2)+L_6*(Beta_1_m2(i-1,2)-Est_Beta_1(i-1));
end

r6(T_ini:end)=Beta_1_m2(T_ini:end,2)-Est_Beta_1(T_ini:end);

%% Residual 7

r7(T_ini:end)=Beta_2_m1(T_ini:end,2)-Beta_2_m2(T_ini:end,2);

%% Residuo 8
Est_Beta_2(1:2,1)=Beta_2_m2(1:2,2);

for i=3:N,
    Est_Beta_2(i)=a8_1*Est_Beta_2(i-1)+a8_2*Est_Beta_2(i-2)+b8_1*Beta_r(i-1,2)+b8_2*Beta_r(i-2,2)+L_8*(Beta_2_m2(i-1,2)-Est_Beta_2(i-1));
end

r8(T_ini:end)=Beta_2_m2(T_ini:end,2)-Est_Beta_2(T_ini:end);

%% Residu 9
r9(T_ini:end)=Beta_3_m1(T_ini:end,2)-Beta_3_m2(T_ini:end,2);

%% Residual 10
Est_Beta_3(1:2,1)=Beta_3_m2(1:2,2);

for i=3:N,
    Est_Beta_3(i)=a10_1*Est_Beta_3(i-1)+a10_2*Est_Beta_3(i-2)+b10_1*Beta_r(i-1,2)+b10_2*Beta_r(i-2,2)+L_10*(Beta_3_m2(i-1,2)-Est_Beta_3(i-1));
end

r10(T_ini:end)=Beta_3_m2(T_ini:end,2)-Est_Beta_3(T_ini:end);

%% Residual 11

Est_Tau_g(1)=Tau_g_m(1,2);

for i=2:N,
    Est_Tau_g(i)=a11_1*Est_Tau_g(i-1)+b11_1*T_g_r(i-1,2)+L_11*(Tau_g_m(i-1,2)-Est_Tau_g(i-1));
end

r11(T_ini:end)=Tau_g_m(T_ini:end,2)-Est_Tau_g(T_ini:end);


%% Residual 12

eta=0.98;
Est_P_g(:,1)=eta*Omega_g_m1(:,2).*Tau_g_m(:,2);

r12(T_ini:end)=P_g_m(T_ini:end,2)-Est_P_g(T_ini:end);

%% Fault detection

for i=1:N,
    
    if abs(r1(i))>th_1,
        fi1(i)=1;
    end
    
    if abs(r2(i))>th_2,
        fi2(i)=1;
    end
    
    if abs(r3(i))>th_3,
        fi3(i)=1;
    end
    
    if abs(r4(i))>th_4,
        fi4(i)=1;
    end
    
    if abs(r5(i))>th_5,
        fi5(i)=1;
    end

    if abs(r6(i))>th_6,
        fi6(i)=1;
    end

    if abs(r7(i))>th_7,
        fi7(i)=1;
    end
    
    if abs(r8(i))>th_8,
        fi8(i)=1;
    end
    if abs(r9(i))>th_9,
        fi9(i)=1;
    end
    
    if abs(r10(i))>th_10,
        fi10(i)=1;
    end
    if abs(r11(i))>th_11,
        fi11(i)=1;
    end
    if abs(r12(i))>th_12,
        fi12(i)=1;
    end
    
end


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




% Fault signals

figure (4)

subplot (4,3,1)
hold on
plot (t,fi1,'b')
title ('fi_1')

subplot (4,3,2)
hold on
plot (t,fi2,'b')
title ('fi_2')

subplot (4,3,3)
hold on
plot (t,fi3,'b')
title ('fi_3')

subplot (4,3,4)
hold on
plot (t,fi4,'b')
title ('fi_4')

subplot (4,3,5)
hold on
plot (t,fi5,'b')
title ('fi_5')

subplot (4,3,6)
hold on
plot (t,fi6,'b')
title ('fi_6')

subplot (4,3,7)
hold on
plot (t,fi7,'b')
title ('fi_7')

subplot (4,3,8)
hold on
plot (t,fi8,'b')
title ('fi_8')

subplot (4,3,9)
hold on
plot (t,fi9,'b')
title ('fi_9')

subplot (4,3,10)
hold on
plot (t,fi10,'b')
title ('fi_{10}')

subplot (4,3,11)
hold on
plot (t,fi11,'b')
title ('fi_{11}')

subplot (4,3,12)
hold on
plot (t,fi12,'b')
title ('fi_{12}')

%% FSM matrix directly taken from the assignment with x and blank spaces
FSM = [
    0 0 0 1 1 0 0 0;
    1 1 1 0 1 1 1 1;
    0 0 0 0 1 0 0 0;
    1 1 1 0 1 1 1 1;
    1 0 0 0 0 0 0 0;
    1 0 0 0 0 0 0 0;
    0 1 0 0 0 0 0 0;
    0 1 0 0 0 1 0 0;
    0 0 1 0 0 0 0 0;
    0 0 1 0 0 0 1 0;
    0 0 0 0 0 0 0 1;
    0 0 0 0 1 0 0 0
    ];

%Observed faults from fig.4
Phi = zeros(8,12);

Phi(1,[5 6]) = 1;
Phi(2,[7 8]) = 1;
Phi(3,[9 10]) = 1;
Phi(4,1) = 1;
Phi(5,[1 2 3 4]) = 1;
Phi(6,8) = 1;
Phi(7,10) = 1;
Phi(8,11) = 1;

% Hamming distance calculation
H = zeros(8,8);

for i = 1:8
    for j = 1:8
        H(i,j) = sum(abs(Phi(i,:) - FSM(:,j)'));
    end
end

disp('Hamming Distance Matrix:')
disp(H)

% Fault isolation 
for i = 1:8

    d_min = min(H(i,:));
    fault = find(H(i,:) == d_min);

    fprintf('Fault %d: Minimum distance = %d, Possible Faults = ',i,d_min);
    disp(fault)

end

%% FSM01

zvf = ones(8,8);
factor01 = zeros(8,8);

for i = 1:8
    for j = 1:8
        for k = 1:12
            if Phi(i,k) == 1 && FSM(k,j) == 0
                zvf(i,j) = 0;
            end
        end

        factor01(i,j) = sum(Phi(i,:).*FSM(:,j)') / sum(FSM(:,j)) * zvf(i,j);

    end
end

disp('zvf Matrix:')
disp(zvf)

disp('FSM01 Factors:')
disp(factor01)

% Isolation Results

for i = 1:8

    factor_max = max(factor01(i,:));
    faults = find(factor01(i,:) == factor_max);

    fprintf('Fault %d: Maximum factor = %.4f, Most probable fault = ',i,factor_max);
    disp(faults)

end