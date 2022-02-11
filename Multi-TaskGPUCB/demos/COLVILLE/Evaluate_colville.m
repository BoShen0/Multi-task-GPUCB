function [Y_eval]= Evaluate_colville(nextPt, nextTasksBatch)

Y_all = -[colville(nextPt,[100 90 10.1 19.8]); 
          colville(nextPt,[90 80 10.1 19.8]);
          colville(nextPt,[100 90 9.1 18]);
          colville(nextPt,[90 80 9.1 18])];
Y_eval = Y_all(nextTasksBatch,:);
end