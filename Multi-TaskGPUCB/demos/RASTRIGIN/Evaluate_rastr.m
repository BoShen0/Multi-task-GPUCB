function [Y_eval]= Evaluate_rastr(nextPt, nextTasksBatch)

Y_all = -[rastr(nextPt,[8 10]); rastr(nextPt,[10 8]);rastr(nextPt,[9 9])];
Y_eval = Y_all(nextTasksBatch,:);
end