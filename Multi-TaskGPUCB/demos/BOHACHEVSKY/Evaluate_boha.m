function [Y_eval]= Evaluate_boha(nextPt, nextTasksBatch)

Y_all = -[boha1(nextPt); boha2(nextPt); boha3(nextPt)];
Y_eval = Y_all(nextTasksBatch,:);
end