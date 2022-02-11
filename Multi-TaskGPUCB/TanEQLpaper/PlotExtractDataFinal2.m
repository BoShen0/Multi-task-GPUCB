
%Shubert Function + Schwefel Function
noiselevel=0.1; RUNS='10-30'; F=100; 
functionname='Shubert Function + Schwefel Function';
functionname2='Shubert Function , Schwefel Function';
Filename='ComparisonSimulationTrial'; 
if(noiselevel==0.1)
    label31='(a),'; label32='(b),'; 
elseif(noiselevel==0.15)
    label31='(c),'; label32='(d),';
end

for run=1:2
if(run==1)
    load(strcat(Filename,num2str(1),'noise',num2str(noiselevel),'.mat'),'x1','x2','trueEQL')
    figure(100+run),mesh(x1,x2,reshape(trueEQL,21,21)');
    xlabel('x_1'),ylabel('x_2'),zlabel('EQL'),xticks([0 0.25 0.5 0.75 1]),yticks([0 0.25 0.5 0.75 1])    
    title(strcat(label31,' y=[',functionname2,']',', target=T_',num2str(run),', noise level=',num2str(noiselevel)))
elseif(run==2)
    load(strcat(Filename,num2str(21),'noise',num2str(noiselevel),'.mat'),'x1','x2','trueEQL')   
    figure(100+run),mesh(x1,x2,reshape(trueEQL,21,21)');
    xlabel('x_1'),ylabel('x_2'),zlabel('EQL'),xticks([0 0.25 0.5 0.75 1]),yticks([0 0.25 0.5 0.75 1])    
    title(strcat(label32,' y=[',functionname2,']',', target=T_',num2str(run),', noise level=',num2str(noiselevel)))  
end
FigureName=strcat(num2str(100+run),'(',RUNS,')',{', '},num2str(noiselevel),' noise.fig');
savefig(string(FigureName))
end

load(strcat(Filename,num2str(40),'noise',num2str(noiselevel),'.mat'))
Per11=zeros(nadd,4,Nrun); Per12=zeros(nadd,4,Nrun); Per13=zeros(nadd,4,Nrun); Per14=zeros(nadd,4,Nrun); Per2=zeros(1,4,Nrun);
for i=1:Nrun
    Per11(:,:,i)=Performance11{i}; Per12(:,:,i)=Performance12{i}; Per13(:,:,i)=Performance13{i}; Per14(:,:,i)=Performance14{i}; Per2(:,:,i)=Performance2{i};  
end
for run=1:2
MP11end=mean(Per11(end,2,((run-1)*nodesigns+1):(run*nodesigns)),3); SP11end=std(Per11(end,2,((run-1)*nodesigns+1):(run*nodesigns)),0,3);
MP12end=mean(Per12(end,2,((run-1)*nodesigns+1):(run*nodesigns)),3); SP12end=std(Per12(end,2,((run-1)*nodesigns+1):(run*nodesigns)),0,3);
MP13end=mean(Per13(end,2,((run-1)*nodesigns+1):(run*nodesigns)),3); SP13end=std(Per13(end,2,((run-1)*nodesigns+1):(run*nodesigns)),0,3);
MP14end=mean(Per14(end,2,((run-1)*nodesigns+1):(run*nodesigns)),3); SP14end=std(Per14(end,2,((run-1)*nodesigns+1):(run*nodesigns)),0,3);
MP2=mean(Per2(end,2,((run-1)*nodesigns+1):(run*nodesigns)),3); SP2=std(Per2(end,2,((run-1)*nodesigns+1):(run*nodesigns)),0,3);
[MP11end SP11end; MP12end SP12end; MP13end SP13end; MP14end SP14end; MP2 SP2]
end

load(strcat('IRSimulationTrial40noise',num2str(noiselevel),'.mat'))
Per15=zeros(nadd,4,Nrun); 
for i=1:Nrun
    Per15(:,:,i)=Performance15{i};  
end

for run=1:2

MP11=mean(Per11(:,2,((run-1)*nodesigns+1):(run*nodesigns)),3); Diff11=Per15(:,2,((run-1)*nodesigns+1):(run*nodesigns))-Per11(:,2,((run-1)*nodesigns+1):(run*nodesigns));
MP12=mean(Per12(:,2,((run-1)*nodesigns+1):(run*nodesigns)),3); Diff12=Per15(:,2,((run-1)*nodesigns+1):(run*nodesigns))-Per12(:,2,((run-1)*nodesigns+1):(run*nodesigns));
MP13=mean(Per13(:,2,((run-1)*nodesigns+1):(run*nodesigns)),3); Diff13=Per15(:,2,((run-1)*nodesigns+1):(run*nodesigns))-Per13(:,2,((run-1)*nodesigns+1):(run*nodesigns));
MP14=mean(Per14(:,2,((run-1)*nodesigns+1):(run*nodesigns)),3); Diff14=Per15(:,2,((run-1)*nodesigns+1):(run*nodesigns))-Per14(:,2,((run-1)*nodesigns+1):(run*nodesigns));
MP15=mean(Per15(:,2,((run-1)*nodesigns+1):(run*nodesigns)),3);
MP2=mean(Per2(end,2,((run-1)*nodesigns+1):(run*nodesigns)),3);

NPDiff11=sum(Diff11>0,3); 
NPDiff12=sum(Diff12>0,3); 
NPDiff13=sum(Diff13>0,3); 
NPDiff14=sum(Diff14>0,3); 

figure(1),hold on,
if(run==1)
    subplot1=subplot(2,2,2*(run-1)+1),title(strcat('(a)',', target=T_',num2str(run))),hold on,    
else
    subplot1=subplot(2,2,2*(run-1)+1),title(strcat('(c)',', target=T_',num2str(run))),hold on,  
end

subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP11,'x','color','blue')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP12,'o','color',[0 0.6 0.2])
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP13,'^','color','red')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP14,'*','color','black')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP15,'square','color','magenta')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP11,'color','blue')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP12,'color',[0 0.6 0.2])
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP13,'color','red')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP14,'color','black')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MP15,'color','magenta')
subplot(2,2,2*(run-1)+1),hold on,plot([1 nadd],[MP2 MP2],':','color','cyan','LineWidth',2)
MIN=round(min([MP11;MP12;MP13;MP14;MP15;MP2])*F)/F;
MAX=round(max([MP11;MP12;MP13;MP14;MP15;MP2])*F)/F;
ylim([MIN 1.4*MAX])
xlabel('Number of added points'), grid on
ylabel('Sample mean of L_1','FontSize',10)   
legend1=legend('EIQ','EMIQ','Q_{0.025}','EIBPQ','EI_I');
set(legend1,'Location','northeast','Orientation','horizontal','NumColumns',2,'FontSize',6,'box','off','FontWeight','bold');

figure(1),hold on,
if(run==1)
    subplot2=subplot(2,2,2*run),title(strcat('(b)',', target=T_',num2str(run))),hold on,    
else
    subplot2=subplot(2,2,2*run),title(strcat('(d)',', target=T_',num2str(run))),hold on,    
end

subplot(2,2,2*run),hold on,plot(1:nadd,NPDiff11,'x','color','blue')
subplot(2,2,2*run),hold on,plot(1:nadd,NPDiff12,'o','color',[0 0.6 0.2])
subplot(2,2,2*run),hold on,plot(1:nadd,NPDiff13,'^','color','red')
subplot(2,2,2*run),hold on,plot(1:nadd,NPDiff14,'*','color','black')
subplot(2,2,2*run),hold on,plot(1:nadd,NPDiff11,'color','blue')
subplot(2,2,2*run),hold on,plot(1:nadd,NPDiff12,'color',[0 0.6 0.2])
subplot(2,2,2*run),hold on,plot(1:nadd,NPDiff13,'color','red')
subplot(2,2,2*run),hold on,plot(1:nadd,NPDiff14,'color','black')
MIN=min([NPDiff11;NPDiff12;NPDiff13;NPDiff14]);
if(mod(MIN,2)==1); MIN=MIN-1; end 
MAX=max([NPDiff11;NPDiff12;NPDiff13;NPDiff14]);
if(mod(MAX,2)==1); MAX=MAX+1; end
yticks(subplot2,[MIN:2:MAX]);
xlabel('Number of added points'),grid on
ylabel('NP_1','FontSize',10)  
legend1=legend('EIQ','EMIQ','Q_{0.025}','EIBPQ');
set(legend1,'Location','southeast','Orientation','horizontal','NumColumns',2,'FontSize',6,'box','off','FontWeight','bold');
if(run==2)
FigureName=strcat(num2str(1),'(',RUNS,')',{', '},num2str(noiselevel),' noise.fig');
savefig(string(FigureName))
end

MXP11=mean(Per11(:,4,((run-1)*nodesigns+1):(run*nodesigns)),3); XDiff11=Per15(:,4,((run-1)*nodesigns+1):(run*nodesigns))-Per11(:,4,((run-1)*nodesigns+1):(run*nodesigns));
MXP12=mean(Per12(:,4,((run-1)*nodesigns+1):(run*nodesigns)),3); XDiff12=Per15(:,4,((run-1)*nodesigns+1):(run*nodesigns))-Per12(:,4,((run-1)*nodesigns+1):(run*nodesigns));
MXP13=mean(Per13(:,4,((run-1)*nodesigns+1):(run*nodesigns)),3); XDiff13=Per15(:,4,((run-1)*nodesigns+1):(run*nodesigns))-Per13(:,4,((run-1)*nodesigns+1):(run*nodesigns));
MXP14=mean(Per14(:,4,((run-1)*nodesigns+1):(run*nodesigns)),3); XDiff14=Per15(:,4,((run-1)*nodesigns+1):(run*nodesigns))-Per14(:,4,((run-1)*nodesigns+1):(run*nodesigns));
MXP15=mean(Per15(:,4,((run-1)*nodesigns+1):(run*nodesigns)),3);

NPXDiff11=sum(XDiff11>0,3); 
NPXDiff12=sum(XDiff12>0,3); 
NPXDiff13=sum(XDiff13>0,3);
NPXDiff14=sum(XDiff14>0,3); 

figure(5),hold on,
if(run==1)
    subplot5=subplot(2,2,2*(run-1)+1),title(strcat('(a)',', target=T_',num2str(run))),hold on,    
else
    subplot5=subplot(2,2,2*(run-1)+1),title(strcat('(c)',', target=T_',num2str(run))),hold on,    
end
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP11,'x','color','blue')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP12,'o','color',[0 0.6 0.2])
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP13,'^','color','red')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP14,'*','color','black')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP15,'square','color','magenta')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP11,'color','blue')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP12,'color',[0 0.6 0.2])
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP13,'color','red')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP14,'color','black')
subplot(2,2,2*(run-1)+1),hold on,plot(1:nadd,MXP15,'color','magenta')
MIN=min([MXP11;MXP12;MXP13;MXP14;MXP15]);
MAX=max([MXP11;MXP12;MXP13;MXP14;MXP15]);
ylim([MIN 1.4*MAX])
legend2=legend('EIQ','EMIQ','Q_{0.025}','EIBPQ','EI_I');
set(legend2,'Location','northeast','Orientation','horizontal','NumColumns',2,'FontSize',6,'box','off','FontWeight','bold');
xlabel('Number of added points'), grid on
ylabel('Sample mean of L_2','FontSize',10)

figure(5),hold on,
if(run==1)
    subplot6=subplot(2,2,2*run),title(strcat('(b)',', target=T_',num2str(run))),hold on,    
else
    subplot6=subplot(2,2,2*run),title(strcat('(d)',', target=T_',num2str(run))),hold on,   
end
subplot(2,2,2*run),hold on,plot(1:nadd,NPXDiff11,'x','color','blue')
subplot(2,2,2*run),hold on,plot(1:nadd,NPXDiff12,'o','color',[0 0.6 0.2])
subplot(2,2,2*run),hold on,plot(1:nadd,NPXDiff13,'^','color','red')
subplot(2,2,2*run),hold on,plot(1:nadd,NPXDiff14,'*','color','black')
subplot(2,2,2*run),hold on,plot(1:nadd,NPXDiff11,'color','blue')
subplot(2,2,2*run),hold on,plot(1:nadd,NPXDiff12,'color',[0 0.6 0.2])
subplot(2,2,2*run),hold on,plot(1:nadd,NPXDiff13,'color','red')
subplot(2,2,2*run),hold on,plot(1:nadd,NPXDiff14,'color','black')
MIN=min([NPXDiff11;NPXDiff12;NPXDiff13;NPXDiff14]); 
if(mod(MIN,2)==1); MIN=MIN-1; end 
MAX=max([NPXDiff11;NPXDiff12;NPXDiff13;NPXDiff14]);
if(mod(MAX,2)==1); MAX=MAX+1; end
yticks(subplot6,[MIN:2:MAX]);
xlabel('Number of added points'),grid on
legend2=legend('EIQ','EMIQ','Q_{0.025}','EIBPQ');
set(legend2,'Location','southeast','Orientation','horizontal','NumColumns',2,'FontSize',6,'box','off','FontWeight','bold');
ylabel('NP_2','FontSize',10)
if(run==2)
FigureName=strcat(num2str(5),'(',RUNS,')',{', '},num2str(noiselevel),' noise.fig');
savefig(string(FigureName))
end

end
