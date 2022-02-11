function [Y_eval]= Evaluate_powell(nextPt, nextTasksBatch)

Y_all = -[powell(nextPt,[10 5 2 10]); 
          powell(nextPt,[10 5 3 9]);
          powell(nextPt,[8 6 2 10])];
Y_eval = Y_all(nextTasksBatch,:);
end