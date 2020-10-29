%% Code for Homework 2, Two-layer perceptron(2020)
% Author: Arshad Nowsath
clc, clear all

%% Loading the given training and validation data from OpenTA
training_data = csvread('training_set.csv');
validation_data = csvread('validation_set.csv');

%% Splitting the data for training and validation set
training_x = training_data(:,1:2);   % contains input patterns 1&2 columns
training_y = training_data(:,3);     % contains targets 3 columns

validation_x = validation_data(:,1:2);
validation_y = validation_data(:,3);

%% Parameters
M1 = 10;     % M1 and M2 are given by random experimental values
M2 = 4;
learning_rate = 0.01; 
epoch = 1000;
classification_error_criteria = 0.12;   % C is below 12%

%% Intialize the weight with some Random Gaussian number
W1 = normrnd(0,1,[M1,2]);
W2 = normrnd(0,1,[M2,M1]);
W3 = normrnd(0,1,[1,M2]);

%% Setting the thresholds to zero
theta_1 = zeros(1,M1);
theta_2 = zeros(1,M2);
theta_3 = zeros(1,1);

%% loop for training the perceptron

for i = 1:epoch
    for j = 1:length(training_x)
        
        mu = randi([1,length(training_y)]);
        % forward propagation
        V1 = tanh(-theta_1' + W1*training_x(mu,:)');
        V2 = tanh(-theta_2' + W2*V1);
        Output = tanh(-theta_3 + W3*V2);
        
        % Backward propagation
        delta_3 = (training_y(mu)-Output)*(1-Output.^2);
        delta_2 = delta_3 * W3.*(1-V2.^2)';
        delta_1 = delta_2*W2.*(1-V1.^2)';
        
        % weights Update
        W3 = W3 + learning_rate * delta_3*V2';
        W2 = W2 + learning_rate * (V1*delta_2)';
        W1 = W1 + learning_rate * delta_1'.*training_x(mu,:);
        
        % Bias Update
        theta_1 = theta_1 - learning_rate * delta_1;
        theta_2 = theta_2 - learning_rate * delta_2;
        theta_3 = theta_3 - learning_rate * delta_3;
    end
    
    % Compute the classification error from the validation set
    V1_validation = tanh(-theta_1 + (W1*validation_x')');
    V2_validation = tanh(-theta_2 + (W2*V1_validation')');
    Output_validation = tanh(-theta_3 + (W3*V2_validation')');
    
    % classification error 
    C = (1/(2*length(validation_y)))*sum(abs(sign(Output_validation)-validation_y));
    
    % Break if the classification error is met
    if C < classification_error_criteria
        disp(['Epoch: ',num2str(i),' C: ',num2str(C)])
        disp('Classification error criteria is met');
        break
    else
        disp(['Epoch: ',num2str(i),' C: ',num2str(C)])
    end     
end

%% Export  weight matrices(W1,W2), vectors(W3) and thresholds(theta) to CSV format

% Weight Matrices
csvwrite('w1.csv',W1);
csvwrite('w2.csv',W2);

% Weight Vectors
csvwrite('w3.csv',W3);

%Thresholds
csvwrite('t1.csv',theta_1);
csvwrite('t2.csv',theta_2);
csvwrite('t3.csv',theta_3);