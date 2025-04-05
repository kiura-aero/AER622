% AER 622 
% Assignment 4
% Kiura Kirubakaran
% 501165598

%% Question 1
x = linspace(0, 1, 100); % chord-wise position

% pressure coefficient equations
Cp_upper = 1 - 15*x.^0.75 .* (1 - x.^0.2) - 0.4*x;
Cp_lower = 1 - 6.4687*x.^0.75 .* (1 - x.^0.2) - 0.4*x;

% pot the pressure coefficients
figure;
plot(x, Cp_upper, 'color', '#4DBEEE','LineWidth', 2); hold on;
plot(x, Cp_lower, 'color', '#0072BD', 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('C_p');
legend('Upper Surface', 'Lower Surface');
title('Incompressible Pressure Coefficient Distribution');
hold off;

%% Question 2 
% Plot the shape of the airfoil that produces these pressure coefficients. 
% The shape of the airfoil that produces these pressure coefficients, and 
% the angle of attack that it is flying at, can be described as follows:

k = 0.1027; % scaling factor

% air foil camber line
y_upper = (k/2).*(15*x.^0.75 .* (1 - x.^0.2));
y_lower = (-k/2).*(6.4687*x.^0.75 .* (1 - x.^0.2));

AOA = 2.86; % angle of attack
theta = deg2rad(AOA); 

% compute rotated coordinates
x_rotated = x * cos(theta) - y_upper * sin(theta);
y_upper_rotated = x * sin(theta) + y_upper * cos(theta);

x_rotated_lower = x * cos(theta) - y_lower * sin(theta);
y_lower_rotated = x * sin(theta) + y_lower * cos(theta);

% plot
figure;
plot(x_rotated, y_upper_rotated, 'color', '#4DBEEE', 'LineWidth', 2); hold on;
plot(x_rotated_lower, y_lower_rotated, 'color', '#0072BD', 'LineWidth', 2);

grid on;
xlabel('x');
ylabel('y');
legend('Upper Surface', 'Lower Surface');
title(['Airfoil at AoA = ', num2str(AOA), 'degrees']);
axis equal;
hold off;


%% Question 3 
% Assuming that we will be cruising at a COMPRESSIBLE Mach number of M=0.80, we 
% apply the Prandtl-Glauert correction to the INCOMPRESSIBLE pressure coefficients, plot
% the COMPRESSIBLE pressure coefficients on the upper and lower surfaces for Mach
% 0.80 flow, and then calculate the Incompressible lift coefficient and the Prandtl-Glauert
% corrected lift coefficient for the Mach 0.80 flow.

M_cruise = 0.80; % COMPRESSIBLE Mach number

% prandtl-gladert correction is being divided into the previous variables
PG_cor = sqrt(1 - M_cruise^2);
Cp_upper_comp = Cp_upper./PG_cor;
Cp_lower_comp = Cp_lower./PG_cor;

% plot
figure;
plot(x, Cp_upper_comp, 'color', '#4DBEEE','LineWidth', 2); hold on;
plot(x, Cp_lower_comp, 'color', '#0072BD', 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('C_p');
legend('Upper Surface', 'Lower Surface');
title('Compressible Pressure Coefficient Distribution');
hold off;

% calculate incompressible lift coeficient
% same cp but made integratable
Cp_upper = @(x) 1 - 15*x.^0.75 .* (1 - x.^0.2) - 0.4*x;
Cp_lower = @(x) 1 - 6.4687*x.^0.75 .* (1 - x.^0.2) - 0.4*x;

% integrand (Cp_upper - Cp_lower)
integrand = @(x) Cp_upper(x) - Cp_lower(x);

% numerical integration to compute the lift coefficient (from lec 14)
C_L = -1 * integral(integrand, 0, 1);  % integral function for numerical integration

disp(['Incompressible Lift Coefficient: ', num2str(C_L)]);

% calculate compressible lift coeficient
% taking cl & dividing by she prandtl on my gladert till im compressible
C_L_comp = C_L ./ PG_cor;

disp(['Compressible Lift Coefficient: ', num2str(C_L_comp)]);

%% Question 4 
% Calculate the critical pressure coefficient for the Mach 0.80 flow and 
% then plot this limit on top of your compressible pressure coefficients. 
% Show that a shock occurs on your upper surface.

% critical pressure coefficient is point where shocks occur. first find
% cp_cr, all point that hit this value have shocks

Y = 1.4; % gamma for air
% equation from lecture 14 (says critical mach number but its for cp)
Cp_critical = (2 /(Y * M_cruise^2))*([(1 + ((Y - 1)/2) * M_cruise^2)/(1 + (Y - 1)/2)]^(Y/(Y-1)) - 1);

disp(['Critical Pressure Coefficient: ', num2str(Cp_critical)]);

% plot
figure;
plot(x, Cp_critical*ones(size(x)), 'r', 'LineWidth', 2, 'LineStyle', '--');hold on;
plot(x, Cp_upper_comp, 'color', '#4DBEEE','LineWidth', 2); 
plot(x, Cp_lower_comp, 'color', '#0072BD', 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('C_p');
legend('Critical Pressure Coefficient', 'Upper Surface', 'Lower Surface');
title('Compressible Pressure Coefficient Distribution with Critical Pressure Coefficient', 'FontName');
hold off;

%% Question 5 
% Plot the critical pressure coefficient for the Mach 0.80 flow on your 
% incompressible pressure coefficients to show that your design’s Sonic 
% Plateau satisfies NASA’s Supercritical Design Guidelines to create an 
% airfoil that is optimized from incompressible to transonic flows at Mach 
% 0.8.

% i jus put again cause matlab wont listen to me
x = linspace(0, 1, 100); % chord-wise position

% pressure coefficient equations
Cp_upper = 1 - 15*x.^0.75 .* (1 - x.^0.2) - 0.4*x;
Cp_lower = 1 - 6.4687*x.^0.75 .* (1 - x.^0.2) - 0.4*x;

% put cp_critical on incomp now
figure;
plot(x, Cp_critical*ones(size(x)), 'r', 'LineWidth', 2, 'LineStyle', '--');hold on;
plot(x, Cp_upper, 'color', '#4DBEEE','LineWidth', 2); 
plot(x, Cp_lower, 'color', '#0072BD', 'LineWidth', 2);
grid on;
xlabel('x');
ylabel('C_p');
legend('Critical Pressure Coefficient', 'Upper Surface', 'Lower Surface');
title('Incompressible Pressure Coefficient Distribution with Critical Pressure Coefficient');
hold off;

%%%
% *Sonic plateau satisfies NASA's Supercritical Design Guidelines*: When
% looking at the graph, the upper surface just touches the critical
% pressure coefficient. When superimposing the critical pressure
% coefficient on the plot, we're looking for places that the cp values
% exceed the cp_critical, which then would mean it does not align with NASA's
% supercritical airfoil guidelines. In my graphs, it is proved that
% this is a supercritical airfoil. 

%% Question 6
% Find the X location of the minimum Cp value.

[Cp_min, idx_min] = min(Cp_upper); % Find the min value and its index
x_min_Cp = x(idx_min); % Get the corresponding x location

%%% 
% The minimum upper pressure coefficient value is also the point that
% touches the critical presssure coefficient. This point should show where
% a shock (if one was to occur) would happen.

% Display the result
disp(['Minimum Cp occurs at x = ', num2str(x_min_Cp), ' with Cp = ', num2str(Cp_min)]);

% this is extra but i wanted to plot it
y_upper_min_cp = interp1(x, y_upper, x_min_Cp);
y_lower_min_cp = interp1(x, y_lower, x_min_Cp);

% plot
figure;
plot(x, y_upper, 'color', '#4DBEEE', 'LineWidth', 2); hold on;
plot(x, y_lower,  'color', '#0072BD', 'LineWidth', 2);
plot(x_min_Cp, y_upper_min_cp, 'wo', 'MarkerSize', 5, 'MarkerFaceColor', 'r'); % Upper surface point

grid on;
xlabel('x');
ylabel('y');
legend('Upper Surface', 'Lower Surface', 'Upper Surface Point', 'Lower Surface Point');
title('Airfoil with Possible Shock Points');
axis equal; 
hold off;
% end of extra stuff

%% Question 7
% The flow along the upper surface accelerates to a supersonic speed, but 
% there is not enough length along the upper surface to isentropically 
% decelerate the flow back to subsonic speeds. Therefore a normal shock is 
% likely to be generated at the location of the maximum Mach number on the 
% upper surface. This normal shock would then, instantaneously decelerate 
% the flow to subsonic. Calculate M_shock, the Mach number of the normal 
% shock

%%%
% A normal shock will be at the location of the maximum Mach number on the
% upper surface. This is where pressure is mininum.  

[Cp_min, idx_min] = min(Cp_upper); % Cp_min is at the shock
% given Cp_min from earlier
Cp_shock = Cp_min; % Value at which shock starts
Y = 1.4;

%% Question 8
% Calculate the Cp_subsonic of the subsonic flow that occurs immediately after the normal
% shock.
 
% define function to find Mach from Cp
Cp_fun = @(M_cruise) (2/(Y*M_cruise.^2))*(((1 + ((Y-1)/2)*M_cruise.^2)/(1 + (Y - 1)/2)).^(Y/(Y-1)) - 1);

mach_guess = 1.2; % initial guess
M_shock = fzero(@(M) Cp_fun(M) - Cp_shock, mach_guess);

disp(['Approximated Mach Number at Shock: ', num2str(M_shock)]);

% Mach after shock using normal shock relations
M_2 = sqrt(((Y - 1)*M_shock^2 + 2) / (2*Y*M_shock^2 - (Y - 1)));

% Get Cp_subsonic after shock
Cp_sub = (2/(Y*M_2^2))*(((1 + ((Y-1)/2)*M_2^2)/(1 + (Y - 1)/2)).^(Y/(Y-1)) - 1);

disp(['Cp_subsonic after shock: ', num2str(Cp_sub)]);

%% Question 9
% Estimate the loss in the lift coefficient that is caused by the presence 
% of the transonic shock. Draw a vertical line from the min Cp_shock to the 
% Cp_subsonic value at the end of the normal shock. Draw a horizontal line 
% at the Cp value at the end of the normal shock until it intersects the Cp 
% distribution on the downstream side of the upper surface Cp distribution. 
% Integrate this area. This area is equal to the loss in the lift 
% coefficient that is produce by the shock. Plot this result to show your 
% pressure distributions, with and
% without the shock.

% the x location of the minimum pressure coefficient is where the shock is
% located
x_shock = x_min_Cp; 

% loss area
x_tail = linspace(x_shock, 1, 100); % from shock to trailing edge
Cp_tail = 1 - 15*x_tail.^0.75 .* (1 - x_tail.^0.2) - 0.4*x_tail;
Cp_flat = Cp_sub * ones(size(x_tail)); % constant Cp after shock

loss = Cp_flat - Cp_tail;
loss_CL = -1 * trapz(x_tail, loss); % loss is negative area

disp(['Lift Coefficient Loss due to Shock: ', num2str(loss_CL)]);

% Plot Cp distributions
figure;
hold on;
plot(x, Cp_upper, 'color', '#4DBEEE', 'LineWidth', 2); % Original Cp
plot(x_tail, Cp_flat, 'r--', 'LineWidth', 2); % Post-shock Cp
% colours in the plot to visualize the area
fill([x_tail, fliplr(x_tail)], [Cp_tail, fliplr(Cp_flat)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');

set(gca, 'YDir', 'reverse'); % flip axis cause prettier
xlabel('x'); ylabel('Cp');
legend('Original Cp', 'Post-Shock Cp', 'Loss in Lift Area');
title('Upper Surface Cp Distribution and Lift Loss due to Shock');
grid on;

%% Question 10
% Make up for this lift coefficient loss by "adding" an equivalent amount of lift to the lower
% surface. Modify the velocity distribution on the lower surface to restore your lift coefficient
% to its intended unshocked value. 

% choose “B” to be a whole number (1,2,3 etc..) and then find the right value of “A”
B = 1;
a = linspace(1, 10, 900);
CL_needed = C_L_comp;
Cp_upper = @(x) 1 - 15*x.^0.75 .* (1 - x.^0.2) - 0.4*x;
Cp_lower = @(x) 1 - 6.4687*x.^0.75 .* (1 - x.^0.2) - 0.4*x;

min_error = inf;
A_best = a(1);
x = linspace(0, 1, 100); % chordwise position for all computations 0-100

for A = a
    % lower pressure coff
    Cp_l_new = @(x) 1 - 6.4687*x.^0.75 .* (1 - x.^0.2) .* (1 - A*(x.^B)) - 0.4*x;
    % integrate for lift
    integrand_new = @(x) Cp_upper(x) - Cp_l_new(x);
    % new lift coeff
    CL_new = -1 * integral(integrand_new, 0, 1);
    % find error
    er_CL = abs(CL_new - CL_needed);
    if er_CL < min_error
        min_error = er_CL;
        A_best = A;
    end
    if er_CL < 0.001
        break;
    end
end

% display result
disp(['B = ', num2str(B)]);
disp(['Best A = ', num2str(A_best)]);
disp(['Recovered Lift Coefficient (CL_new) = ', num2str(CL_new)]);
disp(['Target Lift Coefficient (CL_needed) = ', num2str(CL_needed)]);
disp(['Percent Difference = ', num2str(min_error * 100), '%']);

% plot new lower Cp and compare
Cp_lower_new_plot = 1 - 6.4687*x.^0.75 .* (1 - x.^0.2) .* (1 - A_best*x.^B) - 0.4*x;

figure;
plot(x, Cp_upper(x), 'color', '#4DBEEE', 'LineWidth', 2); hold on;
plot(x, Cp_lower(x), 'color', '#0072BD', 'LineWidth', 2);
plot(x, Cp_lower_new_plot, 'r--', 'LineWidth', 2);
legend('Upper Cp', 'Original Lower Cp', 'Modified Lower Cp');
title('Pressure Distribution with Modified Lower Surface');
set(gca, 'YDir', 'reverse'); % Flip y-axis (Cp is negative up)
xlabel('x'); ylabel('Cp'); grid on;

% New airfoil shape
y_upper = (k/2) .* (15*x.^0.75 .* (1 - x.^0.2));
y_lower = (-k/2) .* (6.4687*x.^0.75 .* (1 - x.^0.2));
y_l_new = (-k/2) .* (6.4687*x.^0.75 .* (1 - x.^0.2) .* (1 - A_best*x.^B));

% plot
figure;
plot(x, y_upper, 'color', '#4DBEEE', 'LineWidth', 2);
hold on;
plot(x, y_l_new, 'r--', 'LineWidth', 2);
hold on;
plot(x, y_lower,'color', '#0072BD', 'LineWidth', 2);
legend('Upper Surface', 'New Lower Surface', 'Original Lower Surface');

axis equal;
title('New Airfoil Geometry');
xlabel('x'); ylabel('y'); grid on;

AOA_new = 2.2; % angle of attack in degrees
theta = deg2rad(AOA_new); % convert to radians

% compute rotated coordinates 
x_rotated = x * cos(theta) - y_upper * sin(theta);
y_upper_rotated = x * sin(theta) + y_upper * cos(theta);

x_rotated_lower_new = x * cos(theta) - y_l_new * sin(theta);
y_lower_rotated_new = x * sin(theta) + y_l_new * cos(theta);

% plot new airfoil 
figure;
plot(x_rotated, y_upper_rotated, 'color', '#4DBEEE', 'LineWidth', 2); hold on;
plot(x_rotated_lower_new, y_lower_rotated_new, 'color', '#0072BD', 'LineWidth', 2);

grid on;
xlabel('x');
ylabel('y');
legend('Upper Surface', 'Lower Surface');
title(['Airfoil at AoA = ', num2str(AOA_new), 'degrees']);
axis equal;
hold off;

%%%
% *A = 2.0912, B = 1*
%%%
% The transonic supercritical airfoil was created using the ideal B value
% fo 1. The iterative code finally determines the value of A to be aprox.
% 2.09. This gives a lift coefficient that matches the unshock
% comprerssible lift coefficient. The needed lift coefficient was 0.83334
% and the new A and B values gave a lift coefficient of 0.8335. This gives
% a percent error of 0.016058%.