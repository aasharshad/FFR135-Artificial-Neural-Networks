%% Code for Homework 2, Linear separability of 4-dimensional Boolean functions(2020)
% Author: Arshad Nowsath
clc, clear all

%% Loading the given training and validation data from OpenTA
training_data = csvread("input_data_numeric.csv");
training_data = training_data(:,2:end);

%% Targets 
E = [1, 1, 1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1];
F = [-1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1];
B = [-1, 1, -1, -1, -1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1];
A = [-1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1];
C = [-1, -1, 1, 1, 1, 1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1];
D = [1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1];

targets = [E' F' B' A' C' D'];

%% Parameters
learning_rate = 0.02;
iteration = 10^5;
tries = 10;
Output = zeros(size(targets));
linearly_separable = zeros(size(targets,2),1);

%% loop for training the perceptron
for i = 1:size(targets,2)
    t = targets(:,i);
    
    for k = 1:tries
        W = -0.02 + rand(4,1)*(0.02-(-0.02)); % Weight Iitialization
        theta = -1 + rand*(1-(-1));           
        for j = 1:iteration
            % Forward propagation
            b = sum(-theta + training_data*W);
            output =  tanh(0.5*(-theta + training_data*W));
            % Stochastic gradient descent
            mu = randi([1,size(targets,1)]);
            d_w = learning_rate*(t(mu)-output(mu))*(1-tanh(b*0.5).^2)*0.5* training_data(mu,:)';
            d_theta = -learning_rate*(t(mu)-output(mu))*(1-tanh(b*0.5).^2)*0.5;
            W = W + d_w;
            theta = theta + d_theta;
            % IF to check whether the linearly seperable break or not
            if sign(tanh(0.5*(-theta + training_data*W))) == t
                break
            end
        end
        if sign(tanh(0.5*(-theta + training_data*W))) == t
            break
        end  
    end
    
    Output(:,i) = sign(tanh(0.5*(-theta + training_data*W)));
    if isequal(Output(:,i),t)
        linearly_separable(i) = true;
    else
        linearly_separable(i) = false;
    end   
end
%% To print which target is linearly seperable 
linearly_separable