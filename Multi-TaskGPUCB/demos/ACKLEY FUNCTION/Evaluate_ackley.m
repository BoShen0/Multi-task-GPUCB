function [Y_eval]= Evaluate_ackley(nextPt, nextTasksBatch)

Y_all = -[ackley(nextPt); ackley(nextPt,15,0.3)];
Y_eval = Y_all(nextTasksBatch,:);
end