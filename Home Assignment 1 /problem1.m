clc
clear all

patterns = [12, 24, 48, 70, 100, 120]';  %pattern generation
N = 120; %No of bits
trails = 10^5; %trails
p_error_list= []; 

for i = 1 : length(patterns) %loop for different pattern generation
    p_i = patterns(i);   % takaing individual pattern for differnt run
    error_count = 0;   % initialize error count to zero
    
    for j = 1 : trails  %loop for 10^5 different trails
        
        p = 2 * randi([0, 1], [N, p_i]) - 1;    % step to generate random pattern with +1 or -1 
        
        random_p = randi(p_i);  % selecting random pattern from each actual pattern 
        random_n = randi(N);   % selecting random neuron from the no of bits
        
        W_i = 1/N * p(random_n,:) * p';   % store a set of p random patterns
        %W_i(random_n) = 0;                % setting diagonal weights to zero     
   
        S0 = sign(p(random_n,random_p));      
        S1 = sign(W_i * p(:,random_p));  % Update single randomly chosen neuron
        
        if S1 == 0
            S1 = 1;          % keeping Signum(0) to Signum(1)
        end
        
        % Check dynamics and see if correct
        if S0 ~= S1
            error_count = error_count + 1; 
        end    
    end
    
    p_error = error_count/trails;
    p_error_list = [p_error_list p_error];
end