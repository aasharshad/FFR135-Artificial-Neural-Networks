function C = classification_error(target, outputs)
length_valset = size(target,1);    
C = 1 / length_valset * sum(outputs ~= target);
end