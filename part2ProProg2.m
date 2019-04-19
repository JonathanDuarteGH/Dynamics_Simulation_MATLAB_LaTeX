% Name: Jonathan Duarte
clear all; clc;
%%%% Values of Scalar Knowns (Dimensions)
pi=4.0*atan(1.0); r1=5.625/12; r2=0.75/12; r3=5.625/12; r4=1.3/12; r5=2; %ft
count=1; g_ft=32.174; %ft/s^2
%%%% Computing mass and the rotational inertia about G for all bodies
% Density of Copper %%[lb/ft^3]%%
rho_cu=559.3545/g_ft; rho_steel=490/g_ft;  
% Lengths and Diameters of links
L=2; D1=2.25/12; L1=3/12; D2=0.5/12; L2=2.5/12; D3=2/12; L3=3.8/12; %ft
% Volume of Motor Windings, Worm Motor, and Worm Wheel %%[ft^3]%%
volume_motor=pi*(D1/2)^2*L1; volume_gear=pi*(D2/2)^2*L2;
volume_wheel=pi*(D3/2)^2*L3;
%Mass of Wiper Blade, Motor Windings, Worm Gear, and Worm Wheel %%[slug]%%
weight_wiper_lbf=1.5; m_wb_slug=weight_wiper_lbf/g_ft;
m_mw_slug=volume_motor*rho_cu; m_wg_slug=volume_gear*rho_steel;
m_ww_slug=volume_wheel*rho_steel;
%Computing the rotational inertia about G for all bodies %%[slugs-ft^2]%%
I_gwb=(1/12)*m_wb_slug*L^2; I_gmw=0.5*m_mw_slug*(D1/2)^2; 
I_gwg=0.5*m_wg_slug*(D2/2)^2; I_gww=0.5*m_ww_slug*(D3/2)^2;
%%%% Initializing initial conditions and guesses
t3r=348*pi/180; %-12 guess
t3radi=348*pi/180; %-12 guess
t4r=221*pi/180; %-139 guess
t4radi=221*pi/180; %-139
t5ri=174*pi/180; %-216
t5r=174*pi/180; %-216 guess
    w_max=39.5*(2*pi/60); %rpm to rad/sec
    T_stall=31*(8.851/12); %Nm to ft*lbf
    %w_rpm=0:0.0041364:w_max-0.001; w_rpm_lin=linspace(0,w_max,1000); 
    %w_rpm_lin=0;
    theta2dot=0;
    theta2new = 0;
    %w_rpm_f=theta2d*30/pi;
    t=0; %seconds
    T_load = 5;
    T_driver=T_stall*(1-theta2dot/w_max);
    theta2_intial=1;
%%%% LOOP
for i = 1:360*4
    t2r = i*(pi/180); %t2r=0:0.01:360*pi/180
    [t3r,t4r,t5r]= nmet4matrixmult(r1,r2,r3,r4,t2r,t3r,t4r,t5r,t5ri,t4radi,t3radi);
    theta2(count)=t2r;
    theta3(count)=t3r;
    theta4(count)=t4r;
    theta5(count)=t5r;
    
    %%%%P A R T 1+2%%%%
    [f_g5xp, f_g5yp, h5p, h4p, h3p, f_g5y, f_g5x, h5, h4, h3, r_g5y, r_g5x] =...
    first_second_kc(r5, theta5, r2, r3, r4, t4r, t3r, t2r, count);

    %%%%P A R T 3%%%%
    %%%%
    
    %%%% Defining the Orientation of the Motor Windings Worm Gear and Worm
    %%%% Wheel
    h2=1; 
    
    [Sum_CAP_A5, Sum_CAP_B5] =...
    powerequation_sim(h2, m_wb_slug, f_g5x, theta2, f_g5y, I_gwb, f_g5xp,...
    f_g5yp, h5, h5p, I_gmw, I_gwg, I_gww);
    
    %%%%P A R T 4%%%%
    %%%% Calculation of T_load applied to link 5
    [T_load]=findTload(h5);
    
    %%%% T_load=F_load_i*(h5./abs(h5)); %F_load acting against T_load: 
    %%%% Expecting step size function
    %%%% Calculation of T_driver applied to link 2
    
    %%%% T_load_i=5; %ft*lbf
    %%%% theta2dot=20*(2*pi/60); %rpm to rad/sec
    
    %%%% Calculation of T_driver applied to link 2
    %%%% h5new=-1.1477; Sum_CAP_A5new=3.7219; Sum_CAP_B5new=0.0129;
    [T_driver]=findTdriver(theta2dot,T_stall,w_max);
    T_motor=T_driver;
   
    %%%%P A R T 5%%%%
    %%%% Computing theta2dot for a purely hypothetical case when
    %%%% theta2=60deg and theta2dot=20rpm  %%%Ans 1.47
    %%%% theta2dotdot=(T_motor+T_load_i.*h5new-Sum_CAP_B5new*theta2dot.^2)/Sum_CAP_A5new;
    theta2dotdot=(T_motor+T_load.*h5-Sum_CAP_B5.*theta2dot.^2)./Sum_CAP_A5;
    
    %%%%P A R T 6%%%%
    %%%% Computing theta5dot, deltatheta2, dt, and deltatheta2dot, dsi
    %%%% for plotting theta2dot vs time and theta5dot vs time
    [theta5dot, deltatheta2, dt, deltatheta2dot, dsi] =...
        eulerfcn(theta2_intial, theta2dot, theta2dotdot, h5);

    arraytheta5dot(i) = theta5dot*(30/pi);
    arrayt(i) = t;
    arraytheta2dot(i) = theta2dot*(30/pi);
    arraydt(i) = dt;
    
    t = t+dt;
    theta2dot = theta2dot+deltatheta2dot;
    theta2new = theta2new+deltatheta2;
    %%%% count=count+1;
end
%%%% Convert all thetas back to degrees
t2d=theta2.*(180/pi);
t3d=theta3.*(180/pi);
t4d=theta4.*(180/pi);
t5d=theta5.*(180/pi);

% t=zeros(1,10,360*pi/180);
% dt2d=zeros();
% t2d=zeros();
% t2dd=zeros();
% dt=zeros();
% 
% theta2dotdot=(T_motor+T_load_i.*h5new-Sum_CAP_B5new*theta2dot.^2)/Sum_CAP_A5new;
% t(1)=0;
% dt(1)=0.035;
% t2d(1)=dt(1)*t2dd(1);
% 
% dt2d(i)=t2dd(i-1)*dt(i-1);
% t(i)=t(i-1)+dt(i-1);
% t2d(i)=t2d0+dt2d(i-1);
% Tm(i)=(1-abs(t2d(i))wmax)*Ts);
% t2d0=t2d(i);
% 
% theta2dotdot(i)=(T_motor+T_load_i.*h5-Sum_CAP_B5*theta2dot.^2)/Sum_CAP_A5;
% dt(i)=-t2d0+

% disp('61103: Part 5 Numbers')
% disp('theta2dot rad*sec(-1)'); disp(theta2dot); 
% disp('theta3 (degrees)'); disp(-360+rad2deg(theta3));
% disp('theta4 (degrees)'); disp(-360+rad2deg(theta4));
% disp('theta5 (degrees)'); disp(-360+rad2deg(theta5));
% disp('h3'); disp(h3); disp('h4'); disp(h4); disp('h5'); disp(h5);
% disp('h3p'); disp(h3p); disp('h4p'); disp(h4p); disp('h5p'); disp(h5p);
% disp('A WIPER BLADE slugs ft^{2}'); disp(CAP_A1);
% disp('B WIPER BLADE slugs ft^{2}'); disp(CAP_B1);
% disp('A Motor Windings slugs ft^{2}'); disp(CAP_A2);
% disp('B Motor Windings slugs ft^{2}'); disp(CAP_B2);
% disp('A Worm Gear slugs ft^{2}'); disp(CAP_A3);
% disp('B Worm Gear slugs ft^{2}'); disp(CAP_B3);
% disp('A Worm Wheel slugs ft^{2}'); disp(CAP_A4);
% disp('B Worm Wheel slugs ft^{2}'); disp(CAP_B4);
% disp('Sum of the As slug ft^{2}'); disp(Sum_CAP_A5);
% disp('Sum of the Bs slug ft^{2}'); disp(Sum_CAP_B5);
% disp('Load Torque slug ft^{2}'); disp(T_load_i); disp('Driving Torque'); disp(T_drivern);
% disp('theta 2 double dot rad*sec^(-2)'); disp(theta2dotdot);

%%%% Graphs for Wiper Motor Slection Part 4: Basic Kinematic Analysis
% figure(1)
% plot(t2d,t3d,'r',t2d,t4d,'b',t2d,t5d,'k');
% axis([0 360 160 370])
% xlabel('\theta_2 (Degrees)','FontSize',14);
% ylabel('\theta_3, \theta_4, and \theta_5 (Degrees)','FontSize',14);
% title...
%     ('Windshield Wiper Mechanism for \theta_3, \theta_4, and \theta_5 vs. \theta_2')
% legend({'\theta_3','\theta_4', '\theta_5'},'Location','Southeast','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(2)
% plot(t2d,h3,'r',t2d,h4,'b',t2d,h5,'k');
% axis([0 360 -1.5 1.5])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('h_3, h_4, and h_5 (Dimensionless)','FontSize',14);
% title...
%     ('Windshield Wiper Mechanism for h_3, h_4, and h_5 (Dimensionless) vs. \theta_2')
% legend({'h_3','h_4', 'h_5'},'Location','Southeast','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(3)
% plot(t2d,h3p,'r',t2d,h4p,'b',t2d,h5p,'k');
% axis([0 360 -2 1.5])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('h_3^{\prime}, h_4^{\prime}, and h_5^{\prime} (Dimensionless)','FontSize',14);
% title...
%     ('Windshield Wiper Mechanism for h_3^{\prime},h_4^{\prime}, and h_5^{\prime} vs. \theta_2')
% legend({'h_3^{\prime}','h_4^{\prime}','h_5^{\prime}'},'Location','South','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(4)
% plot(t2d,f_g5x,'r',t2d,f_g5y,'b');
% axis([0 360 -1.5 1.5])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('f_{g5x} and f_{g5y} (Ft)','FontSize',14);
% title...
%     ('First Order KCs for Windshield Wiper Mechanism mass centers vs. \theta_2')
% legend({'f_{g5x}','f_{g5y}'},'Location','South','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(5)
% plot(t2d,f_g5xp,'r',t2d,f_g5yp,'b');
% axis([0 360 -1.5 1.5])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('f_{g5x}^{\prime} and f_{g5y}^{\prime} (Ft)','FontSize',14);
% title...
%     ('Second Order KCs for Windshield Wiper Mechanism mass centers vs. \theta_2')
% legend({'f_{g5x}^{\prime}','f_{g5y}^{\prime}'},'Location','North','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(6)
% plot(t2d,CAP_A1,'r',t2d,CAP_B1,'b');
% axis([0 360 -0.05 0.1])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('A and B (Slugs ft^{2})','FontSize',14);
% title...
%     ('A and B for Wiper Blade vs. \theta_2')
% legend({'A','B'},'Location','North','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(7)
% plot(t2d,CAP_A2,'r',t2d,CAP_B2,'b');
% axis([0 360 0 4])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('A and B (Slugs ft^{2})','FontSize',14);
% title...
%     ('A and B for Motor Windings vs. \theta_2')
% legend({'A','B'},'Location','North','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(8)
% plot(t2d,CAP_A3,'r',t2d,CAP_B3,'b');
% axis([0 360 0 7E-3])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('A and B (Slugs ft^{2})','FontSize',14);
% title...
%     ('A and B for Worm Gear vs. \theta_2')
% legend({'A','B'},'Location','North','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(9)
% plot(t2d,CAP_A4,'r',t2d,CAP_B4,'b');
% axis([0 360 0 4E-4])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('A and B (Slugs ft^{2})','FontSize',14);
% title...
%     ('A and B for Worm Wheel vs. \theta_2')
% legend({'A','B'},'Location','North','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(10)
% plot(t2d,Sum_CAP_A5,'r',t2d,Sum_CAP_B5,'b');
% axis([0 360 -1 5])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('Sum A and Sum B (Slugs ft^{2})','FontSize',14);
% title...
%     ('Sum of As and Sum Bs for Windshield Wiper Mechanism vs. \theta_2')
% legend({'Sum A',' Sum B'},'Location','North','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(11)
% plot(t2d,T_load,'b');
% axis([0 370 -6 6])
% xlabel('\theta_2 (degrees)','FontSize',14);
% ylabel('T_{load} (slug ft^2)','FontSize',14);
% title...
%     ('The Load Torque Acting on Link 5 vs. \theta_2')
% legend({'T_{load}'},'Location','North','FontSize',12);
% set(gca,'FontSize',14);
% 
% figure(12)
% plot(w_rpm_f,T_driver,'b');
% xlabel('$\dot{\theta_2}$(RPM)','interpreter','latex','FontSize',14);
% ylabel('T_{driver} (slug ft^2)','FontSize',14);
% title...
%     ('The Driving Torque Acting on Link 2 vs. \theta_2')
% legend({'T_{driver}'},'Location','North','FontSize',12);
% set(gca,'FontSize',14);

figure(13)
plot(arrayt,arraytheta2dot)
xlabel('time (s)')
ylabel('$$\dot{\theta}_{2}$$ (rpm)','Interpreter','latex')
set(gca, 'FontSize', 16)
 
figure(14)
plot(arrayt,arraytheta5dot)
xlabel('time (s)')
ylabel('$$\dot{\theta}_{5}$$ (rpm)','Interpreter','latex')
set(gca, 'FontSize', 16)
 
theta2dotmax = max(arraytheta2dot);
theta2dotmin = min(arraytheta2dot(1,720:end));
theta2dotavg = (theta2dotmax+theta2dotmin)/2;
C_f = (theta2dotmax-theta2dotmin)/theta2dotavg;
[pks, indexes]=findpeaks(arraytheta5dot);
Tvalues=arrayt(indexes);
Time_req = Tvalues(4)-Tvalues(3);

% disp('61103: Part 5 Numbers')
% disp('At Steady State'); 
% disp('Theta2dot max (rpm)'); disp(theta2dotmax); 
% disp('Theta2dot min (rpm)'); disp(theta2dotmin); 
% disp('Theta2dot Average (rpm)'); disp(theta2dotavg);
% disp('Coefficient of Fluctuation (percentage)'); disp(C_f*100);
% disp('Time required to wipe the windsield clean (seconds)'); disp(Time_req);
% disp('Maxmimum Wiper load the Mechanism can overcome (slug ft^2)'); disp(T_load);