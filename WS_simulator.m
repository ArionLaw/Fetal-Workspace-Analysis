% inputs
% Setup 1
L = [0.5, 0.5, 0];    % [ls, le, lw], in m
rotation = 30;         % base rotation angle
x_shift = 0.8;          % base x position
y_shift = 0.5;          % base y position
Limits_up=[145+rotation 120 0.3];%[theta1 theta2 d1] d2 not needed for x-y
Limits_dw=[0+rotation -120 0.1];%[theta1 theta2 d1] d2 not needed for x-y

%% Configuration Space Sample 
N_rand=10000;

%% % this is a trick I used to generate a grid of values, you can do ir differently
Levels=floor(abs(Limits_up-Limits_dw)/10)+1; 
Offset=(Limits_dw/10-1);
D=fullfact(Levels);
QS=[(D(:,1)+Offset(1))*10 (D(:,2)+Offset(2))*10 (D(:,3)+Offset(3))*10];

%% Here I also enriched the points by adding random points from
% configuration space
QS_rand=[rand(N_rand,1)*(Limits_up(1)-Limits_dw(1))+Limits_dw(1)  rand(N_rand,1)*(Limits_up(2)-Limits_dw(2))+Limits_dw(2) rand(N_rand,1)*(Limits_up(3)-Limits_dw(3))+Limits_dw(3)];
QS=[QS;QS_rand];

%% FK model
s1 = sind(QS(:,1));
s12 = sind(QS(:,1) + QS(:,2));

c1 = cosd(QS(:,1));
c12 = cosd(QS(:,1) + QS(:,2));

% End-effector's positions
x = x_shift + L(1)*c1 + (L(2)+QS(:,3)).*c12;   
y = y_shift + L(1)*s1 + (L(2)+QS(:,3)).*s12;

%% plot the WS
% figure()
% plot(x,y,'.r')
plot(x,y,'c.')
title('Workspace X-Y')
xlabel('x-axis [m]') 
ylabel('y-axis [m]') 
set(gca,'fontsize',14);
hold on
plot(0,0,'+r')

%%
x_tasks = [1.7 -0.25 0];
y_tasks = [1.3 1 0];
labels = {'A','B','C'};
hold on
p = plot(x_tasks,y_tasks,'r*');
p.MarkerSize = 8;
%}