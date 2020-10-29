% Assignment 1 - problem 3
% Stochastic Hopfield network
clear all;close all; clc;

%given parameters 

T = 2*10^5;
N = 200; % neurons to use
p = 7; % random patterns to store

mu = zeros(1,100);   

for i = 1:100     %in order to repeat the experiment for 100 times
   
    rand_patterns = 2 * randi([0, 1], [N, p]) - 1;   % generate random patterns of row N and column p
    
    W = zeros(N, N);  % creating weight matrix of N*N 
    
    for j = 1 : p           %steps to store the patterns in the network
        W = W + rand_patterns(:, j) * rand_patterns(:, j)';  %hebbs rule
    end
        
    W = W / N; % Normalize the weight matrix
    W = W - diag(diag(W)); % diagonal elements to zero, comment to get ans for 3a
    
    S0 = rand_patterns(:,1);
    
    m1_T = zeros(1,T);  %initialise order parameter to zeros of 1*T
    
    for t = 1 : T   %for loop to get the order parameters
        S1 = S0;
        ni = randi(N);   %making random N 
        bi = W(ni,:) * S0;     
        probability = sigmf(bi, [4,0]);
        S1(ni) = randsrc(1, 1, [1,-1; probability, 1-probability]);
        m1_T(t) = 1/N * S1' * rand_patterns(:,1);
        S0 = S1;
    end    
    mu(i) = 1 / T * sum(m1_T);    %order parameter or finite time average
end

m1_T_avg = 1/100 * sum(mu);  %resulting average order parameter
