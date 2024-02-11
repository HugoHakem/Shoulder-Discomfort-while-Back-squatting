clear all

%% Body dimension
H1_L=[0.34 0.37 0.40]; % Forearm (m)
H2_L=[0.28 0.30 0.32]; % Humerus (m)
c_L=[0.025 0.03 0.035]; % Half of the side of the shoulder associated to a cube (m)
bheta_L=[4 10 15]; % Fixed angle in degrees between forearm and humerus in the YX plan
Bheta_L=bheta_L*pi/180; % Fixed angle in Radian between forearm and humerus in the YX plan

H1_opt=H1_L(2);
H2_opt=H2_L(2);
c_opt=c_L(1);
Bheta_opt=Bheta_L(1);

OPTIMUM_DeltaH=zeros(4,3);
OPTIMUM_DeltaHand=zeros(4,3);
OPTIMUM_AnglXY=zeros(4,3);
OPTIMUM_AnglZY=zeros(4,3);
for test=1:1:12
    if test<=3
        row_tested=1;
        col_test=test;
        H1=H1_L(test);
        H2=H2_opt;
        c=c_opt;
        Bheta=Bheta_opt;
    elseif (test>3)&&(test<=6)
        row_tested=2;
        col_test=test-3;
        H1=H1_opt;
        H2=H2_L(test-3);
        c=c_opt;
        Bheta=Bheta_opt;
    elseif (test>6)&&(test<=9)
        row_tested=3;
        col_test=test-6;
        H1=H1_opt;
        H2=H2_opt;
        c=c_L(test-6);
        Bheta=Bheta_opt;
    else
        row_tested=4;
        col_test=test-9;
        H1=H1_opt;
        H2=H2_opt;
        c=c_opt;
        Bheta=Bheta_L(test-9);
    end
%% Moment arm
epsilon=0; % Moment arm of the muscle force in the X direction (m) Not used actually as it accounted for the barbell not being exactly on the shoulders but no need to be tat precise. 
a=0.05; % Moment arm of the muscle force in the YZ plan (m)
b=0.1; % Moment arm of the Reaction force due to the weight in Z direction (m)


%% Variable parameters
n=13; % Number of variable
m=15;
DeltaH=linspace(0.01,2*c,n).'; % Barbell positionment on the shoulder, starting from the top of shoulder to the end, from 0 to 2*c (m)
DeltaHand=linspace(0.10,(H1+H2)*(1-0.1),m); %Spacing of the hand on the barbell, measured between the shoulders middle to the hand, from 2*c (you don't want your hand at the edge of the shoulder to less than H1+H2 (m)

%% Equation solver angle between humerus and forearm
Theta=zeros(n,m);
Gamma=zeros(n,m);
for i=1:1:n
    for j=1:1:m
        try
            syms theta gamma
            eqn1 = H1*cos(gamma)*cos(Bheta)-H2*cos(theta) == c-DeltaH(i,1);
            eqn2 = H1*sin(gamma)*cos(Bheta) + H2*sin(theta) == DeltaHand(1,j);
            eqn3 = theta >= 0 & theta <= pi/2;
            eqn4 =  gamma >= 0 & gamma <= pi/2;
            [theta, gamma]=solve(eqn1,eqn2,eqn3,eqn4,theta,gamma);
            Theta(i,j)=theta;
            Gamma(i,j)=gamma;
        catch 
            Theta(i,j)=nan;
            Gamma(i,j)=nan;
        end
    end
end
%% Test of the manual solving of the equation system but no way to decide on the solutions easily. 
% 
% At=(c-DeltaH.').^2+(H2)^2+DeltaHand.^2-(H1*cos(Bheta))^2-2*H2*(c-DeltaH.');
% Bt=4*H2*DeltaHand;
% Ct=(c-DeltaH.').^2+(H2)^2+DeltaHand.^2-(H1*cos(Bheta))^2+2*H2*(c-DeltaH.');
% Deltat=Bt.^2-4*At.*Ct;
% 
% Theta_bis_p=2*atan((Bt+sqrt(Deltat))./(2*At));
% Theta_bis_m=2*atan((Bt+sqrt(Deltat))./(2*At));
% 
% Ag=(c-DeltaH.').^2+(H1*cos(Bheta))^2+DeltaHand.^2-H2^2+2*H1*cos(Bheta)*(c-DeltaH.');
% Bg=4*H1*cos(Bheta)*DeltaHand;
% Cg=(c-DeltaH.').^2+(H1*cos(Bheta))^2+DeltaHand.^2-H2^2-2*H1*cos(Bheta)*(c-DeltaH.');
% Deltag=Bg.^2-4*Ag.*Cg;
% 
% Gamma_bis_p=2*atan((Bg+sqrt(Deltag))./(2*Ag));
% Gamma_bis_m=2*atan((Bg-sqrt(Deltag))./(2*Ag));
%% Upper Body / shoulder inclination bis
Alpha = f_alpha(DeltaH,Gamma,Theta,c,H1,H2,Bheta);

%% Load on the back due to the barbell 
pho=0.1; %portion of Hand resistance regarding the total reaction force (from 0 to 1)
W=30;
Rhand=f_rhand(pho,W,Alpha,Bheta,Gamma);
Rn=W*cos(Alpha)-Rhand*sin(Bheta).*cos(Gamma);
T=W*sin(Alpha)-Rhand*cos(Bheta).*cos(Gamma);
%% Shoulder forces and Moments

FM=(-Rhand*sin(Bheta).*cos(Gamma).*DeltaHand+Rn*b+Rhand*cos(Bheta).*sin(Gamma)*(c+epsilon))./(sin(Theta)*(c+epsilon));
Fx=-Rn-Rhand*sin(Bheta).*cos(Gamma);
Fy=T+Rhand*cos(Bheta).*cos(Gamma)+FM.*cos(Theta);
Fz=-FM.*sin(Theta)+Rhand*cos(Bheta).*sin(Gamma);
Mx=FM*a+T*b-Rhand*cos(Bheta).*cos(Gamma).*DeltaHand+Rhand*cos(Bheta).*sin(Gamma).*(c-DeltaH);
My=0;
Mz=-T*c-FM.*cos(Theta)*(c+epsilon)+Rn.*(c-DeltaH)+Rhand*sin(Bheta).*cos(Gamma).*(c-DeltaH)-Rhand*cos(Bheta).*cos(Gamma)*c;

%% Plot abs(FM)
% figure;
% subplot(1,2,1)
% im=imagesc(DeltaHand,DeltaH,abs(FM));
% title('Muscle Force FM (N)')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(abs(FM)))
% colorbar
% 
% subplot(1,2,2)
% sortedFM=sort(abs(FM(~isnan(FM))));
% im=imagesc(DeltaHand,DeltaH,abs(FM),[sortedFM(1) sortedFM(10)]);
% title('Muscle Force FM (N)')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(abs(FM)))
% colorbar
%% Plot F
% figure
% im=imagesc(DeltaHand,DeltaH,sqrt(Fx.^2+Fy.^2+Fz.^2));
% title('Joint reaction force (N)')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(abs(FM)))
% colorbar

%% plot Fx , Fy , Fz
% figure
% subplot(1,3,1)
% im=imagesc(DeltaHand,DeltaH,abs(Fx));
% title('Joint reaction force Fx (N)')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(Fx))
% colorbar
% 
% subplot(1,3,2)
% im=imagesc(DeltaHand,DeltaH,abs(Fy));
% title('Joint reaction force Fy (N)')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(Fy))
% colorbar
% 
% subplot(1,3,3)
% im=imagesc(DeltaHand,DeltaH,abs(Fz));
% title('Joint reaction force Fz (N)')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(Fz))
% colorbar

%% Plot Fxy oriantation
% subplot(1,2,1)
% im=imagesc(DeltaHand,DeltaH,sqrt(Fx.^2+Fy.^2));
% title('Joint reaction force Fxy (N)')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(abs(Fx)))
% colorbar
% 
% subplot(1,2,2)
% im=imagesc(DeltaHand,DeltaH,(atan(Fy./Fx)+pi)*180/pi);
% title('Angle between x axis and Fxy in degrees')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(Fx))
% colorbar

%% Plot Fzy oriantation 
% figure
% subplot(1,2,1)
% im=imagesc(DeltaHand,DeltaH,sqrt(Fy.^2+Fz.^2));
% title('Joint reaction force Fzy (N)')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(abs(Fy)))
% colorbar
% 
% subplot(1,2,2)
% im=imagesc(DeltaHand,DeltaH,(atan(Fy./Fz)+pi)*180/pi);
% title('Angle between z axis and Fzy in degrees')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(Fy))
% colorbar

%% Plot Fzx oriantation Nothing really interesting, keep being around 20 degrees
% figure
% subplot(1,2,1)
% im=imagesc(DeltaHand,DeltaH,sqrt(Fz.^2+Fx.^2));
% title('Joint reaction force Fzx (N)')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(abs(Fz)))
% colorbar
% 
% subplot(1,2,2)
% im=imagesc(DeltaHand,DeltaH,(atan(Fx./Fz))*180/pi);
% title('Angle between z axis and Fzx in degrees')
% xlabel('DeltaHand (m)')
% ylabel('DeltaH (m)')
% set(im,'alphadata',~isnan(Fz))
% colorbar

%% Comparison between angle in xy axis and zy axis
ang_xy_max=105;
ang_zy_min=110;

figure
ax1=subplot(1,2,1);
ang_xy=((atan(Fy./Fx)+pi)*180/pi);
im1=imagesc(DeltaHand,DeltaH,ang_xy,[90,ang_xy_max]);
title('Angle between x axis and Fxy in degrees')
subtitle(append('H1 = ',string(H1),' / H2 = ',string(H2),' / c = ',string(c),' / Bheta = ',string(Bheta*180/pi)))
xlabel('DeltaHand (m)')
ylabel('DeltaH (m)')
set(im1,'alphadata',~isnan(Fx))
set(im1,'alphadata',(ang_xy<ang_xy_max)+(ang_xy>ang_xy_max)/4)
colormap(ax1,flip(parula))
colorbar


ax2=subplot(1,2,2);
ang_zy=((atan(Fy./Fz)+pi)*180/pi);
im2=imagesc(DeltaHand,DeltaH,ang_zy,[ang_zy_min,150]);
title('Angle between z axis and Fzy in degrees')
subtitle(append('H1 = ',string(H1),' / H2 = ',string(H2),' / c = ',string(c),' / Bheta = ',string(Bheta*180/pi)))
xlabel('DeltaHand (m)')
ylabel('DeltaH (m)')
set(im2,'alphadata',~isnan(Fy))
set(im2,'alphadata',(ang_zy>ang_zy_min)+(ang_zy<ang_zy_min)/5)
colormap(ax2,'parula')
colorbar

figure;

ang_xy_novar=ang_xy./std(ang_xy,0,2,'omitnan');
ang_zy_novar=ang_zy./std(ang_zy,0,2,'omitnan');
compar_ang_xy_zy=ang_zy_novar./ang_xy_novar;
im3=imagesc(DeltaHand,DeltaH, compar_ang_xy_zy);
title('Comparison Normalized std Angle_{zy}/Angle_{xy}')
subtitle(append('H1 = ',string(H1),' / H2 = ',string(H2),' / c = ',string(c),' / Bheta = ',string(Bheta*180/pi)))
xlabel('DeltaHand (m)')
ylabel('DeltaH (m)')
set(im3,'alphadata',~isnan(Fx))
set(im3,'alphadata',((ang_zy>ang_zy_min)&(ang_xy<ang_xy_max))+((ang_xy>ang_xy_max)|(ang_zy<ang_zy_min))/5)
colorbar

compar_ang_xy_zy(((ang_xy>ang_xy_max)|(ang_zy<ang_zy_min)))=nan;
[i_opt,j_opt]=find(compar_ang_xy_zy==max(max(compar_ang_xy_zy)));
OPTIMUM_DeltaH(row_tested,col_test)=DeltaH(i_opt,1)/(2*c);
OPTIMUM_DeltaHand(row_tested,col_test)=DeltaHand(1,j_opt)/(H1+H2);
OPTIMUM_AnglXY(row_tested,col_test)=ang_xy(i_opt,j_opt);
OPTIMUM_AnglZY(row_tested,col_test)=ang_zy(i_opt,j_opt);

end
%% Box chart 
Parameters=["Forearm"; "Forearm";"Forearm";"Humerus";"Humerus";"Humerus";"Shoulder size";"Shoulder size";"Shoulder size";"Bheta";"Bheta";"Bheta"];
OPTIMUM_DeltaH_line=[];
for i=1:1:4
OPTIMUM_DeltaH_line=[OPTIMUM_DeltaH_line; OPTIMUM_DeltaH(i,:).'];
end
figure
subplot(1,2,1)
boxchart(categorical(Parameters),OPTIMUM_DeltaH_line)
ylabel('Optimal Delta_H / (2*c)')
xlabel('Parameters')

OPTIMUM_DeltaHand_line=[];
for i=1:1:4
OPTIMUM_DeltaHand_line=[OPTIMUM_DeltaHand_line; OPTIMUM_DeltaHand(i,:).'];
end
%figure
subplot(1,2,2)
boxchart(categorical(Parameters),OPTIMUM_DeltaHand_line)
ylabel('Optimal Delta_H_a_n_d / (H1+H2)')
xlabel('Parameters')

OPTIMUM_AnglXY_line=[];
for i=1:1:4
OPTIMUM_AnglXY_line=[OPTIMUM_AnglXY_line; OPTIMUM_AnglXY(i,:).'];
end
figure
subplot(1,2,1)
boxchart(categorical(Parameters),OPTIMUM_AnglXY_line)
ylabel('Optimal AnglXY (°)')
xlabel('Parameters')

OPTIMUM_AnglZY_line=[];
for i=1:1:4
OPTIMUM_AnglZY_line=[OPTIMUM_AnglZY_line; OPTIMUM_AnglZY(i,:).'];
end
%figure
subplot(1,2,2)
boxchart(categorical(Parameters),OPTIMUM_AnglZY_line)
ylabel('Optimal AnglZY (°)')
xlabel('Parameters')