function [Y_eval]= Evaluate_levy(nextPt, nextTasksBatch)

Y_all = -[levy(nextPt,[1 10 1]); levy(nextPt,[1.2 9 1.2])];
Y_eval = Y_all(nextTasksBatch,:);
end