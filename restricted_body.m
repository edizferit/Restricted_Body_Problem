% ME303 TERM PROJECT QUESTION 2
clear all, close all, clc
format long

% 2: RESTRICTED 3-BODY PROBLEM

%================= GIVENS =================%
% MASSES
m1 = 0.012277471;
m2 = 1 - m1;
% PERIOD
T = 17.06521656015796; % Time required for one full round of orbit
% STEP SIZES
h_e = T/24000;
h_rk4 = T/6000;
% INTERVAL
tmin = 0;
tmax = T;
% INITIAL VALUES
yi(1) = 0.994;
yi(2) = 0;
yi(3) = 0;
yi(4) = -2.0015851063790825;

%======== FIRST ORDER SYSTEM OF DIFFERENTIAL EQUATIONS ========%
dy{1} = @(y) y(2);
dy{2} = @(y) y(1) + 2*y(4) - m2*(y(1) + m1)/((y(1) + m1).^2 + y(3).^2)...
    .^(3/2) - m1*(y(1) - m2)/((y(1) - m2).^2 + y(3).^2).^(3/2);
dy{3} = @(y) y(4);
dy{4} = @(y) y(3) - 2*y(2) - m2*y(3)/((y(1) + m1).^2 + y(3).^2).^(3/2)...
    - m1*y(3)/((y(1) - m2).^2 + y(3).^2).^(3/2);


% ======== EULER SOLUTION ========%
% Values For Euler Solution
% STEP SIZE
h = h_e;
% INTERVAL
t = tmin:h:tmax;
% NUMBER OF ITERATIONS
n = (tmax-tmin)/h;
% POSITION VECTORS
y =  zeros(4,n);
% INITIAL CONDITIONS
y(:,1) = yi(:);

for i=1:n
    for j=1:4
        y(j,i+1) = y(j,i) + h*dy{j}(y(:,i));
    end
end
y1 = y(1,:);
y2 = y(3,:);

figure('Position', [100 180 1400 400])
subplot(1,3,1)
plot(y1,y2, 'g', 'linewidth', 1.5)
hold on, grid on
plot(0,0, '.r', 'linewidth', 3)
title('Euler Solution')
xlabel('y1')
ylabel('y2')
legend('Positions', 'Origin');


%======== 4TH ORDER RUNGE KUTTA SOLUTION ========%
% Values For RK4 Solution
% STEP SIZE
h = h_rk4;
% INTERVAL
t = tmin:h:tmax;
% NUMBER OF ITERATIONS
n = (tmax-tmin)/h;
% POSITION VECTORS
y =  zeros(4,n);
% INITIAL CONDITIONS
y(:,1) = yi(:);

f{1} = @(g,y,h,fv) g(y);
f{2} = @(g,y,h,fv) g(y + h/2*fv(:,1));
f{3} = @(g,y,h,fv) g(y + h/2*fv(:,2));
f{4} = @(g,y,h,fv) g(y + h*fv(:,3));
fv = zeros(4,4); % Each row stores f values for each equation
for i=1:n
    for j=1:4 % four f values
        for k=1:4 % four equations
            fv(k,j) = f{j}(dy{k}, y(:,i), h, fv);
        end
    end
    y(:,i+1) = y(:,i) + h/6*(fv(:,1) + 2*fv(:,2) + 2*fv(:,3) + fv(:,4));
end
y1 = y(1,:);
y2 = y(3,:);

subplot(1,3,2)
plot(y1,y2, 'b', 'linewidth', 1.5)
hold on, grid on
plot(0,0, '.r', 'linewidth', 3)
title('RK4 Solution')
xlabel('y1')
ylabel('y2')
legend('Positions', 'Origin');


%=============== ODE45 SOLUTION ===============%
tspan = [tmin tmax];
y0 = [yi(1) yi(2) yi(3) yi(4)];
[t,y] = ode45(@(t,y) odefcn(y,m1,m2), tspan, y0);
y1 = y(:,1);
y2 = y(:,3);

subplot(1,3,3)
plot(y1,y2, 'k', 'linewidth', 1.5)
hold on, grid on
plot(0,0, '.r', 'linewidth', 3)
title('ODE45 Solution')
xlabel('y1')
ylabel('y2')
legend('Positions', 'Origin');

function dydt = odefcn(y,m1,m2)
dydt = zeros(4,1);
dydt(1) = y(2);
dydt(2) = y(1) + 2*y(4) - m2*(y(1) + m1)/((y(1) + m1).^2 + y(3).^2)...
      .^(3/2) - m1*(y(1) - m2)/((y(1) - m2).^2 + y(3).^2).^(3/2);
dydt(3) = y(4);
dydt(4) = y(3) - 2*y(2) - m2*y(3)/((y(1) + m1).^2 + y(3).^2).^(3/2)...
      - m1*y(3)/((y(1) - m2).^2 + y(3).^2).^(3/2);
end

