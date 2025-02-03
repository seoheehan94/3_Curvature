% Define theta values from 0 to pi
theta = [0*pi 0.255*pi 0.452*pi 0.58*pi 0.664*pi 0.722*pi 0.764*pi 0.795*pi];
% Adjust the number of points for smoother plot

% Calculate bA for each theta
% bA_values = 1 ./ tan((pi/2) - (theta_values/2));
bA_values=1./(tan((pi/2)-(theta/2)));
bA_values
% Plot bA as a function of theta
figure;
plot(theta, bA_values, 'bo-', 'LineWidth', 1.5);
xlabel('Theta (radians)');
ylabel('bA');
title('Graph of bA as a function of Theta');
grid on;


xLength=30;
yLength=30;
mA=4;
waveLength=1;
alpha=0;

% Define bA values
%[0, 0.4286, 0.8571, 1.2857, 1.7143, 2.1429, 2.5714, 3.0000];
% theta = [0*pi 0.003*pi 0.007*pi 0.014*pi 0.024*pi 0.045*pi 0.12*pi 0.3*pi];
theta = linspace(0,157.5,8);
theta = theta/1200;
theta(3)=theta(3)*0.9;
theta(4)=theta(4)*0.9;
theta(6)=theta(6)*1.2;
theta(7)=theta(7)*2.4;
theta(8)=theta(8)*5;

% Create subplots
figure;
for i = 1:numel(theta)
    % Generate the filter for current bA
    [SpaceKernel, ~, ~] = AngleFilter(waveLength, alpha, theta(i), mA, xLength, yLength);
    filter = real(SpaceKernel);
    
    % Plot in subplot
    subplot(1, 8, i);
    imagesc(filter);
    axis equal tight;
    title(['theta = ', num2str(theta(i))]);
    colormap('gray');
end
sgtitle('Filters for Different theta Values');