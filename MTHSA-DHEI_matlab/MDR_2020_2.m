% function [err_rate,LR,PSI, CE,ED,F] = MDR_2020_2(SNP_COM,state)
function [CE, CA, PSI,F] = MDR_2020_2(SNP_COM,state)
%% 没有拓**同意，代码不要外传 和 网上发布
%% reference: 
    % William S et al. (2008) Alternative contingency table measures improve the power and detection of multifactor dimensionality reduction 

% SNP_COM 基因型 ： 0 1 2 
% state   表型 1 case ； 0 control

%%
ua = unique(state);
if length(ua)~=2 
    disp(' Class state is not equal to 2 !!!');
elseif min(ua)>0
    state = state -min(ua);
end

[xrow,xcol] = size(SNP_COM);
[Data,idx,cid]=unique(SNP_COM,'rows');
[lrow,~]=size(Data);
sample=zeros(lrow,1);
disease=sample;
control=sample;

for i=1:xrow   %% 统计每个基因型组合出现的次数
   if state(i) == 1
       disease(cid(i)) = disease(cid(i)) + 1;
   else
       control(cid(i)) = control(cid(i)) + 1;
   end
end

ProCaseToControl = (disease+0.1)./(control + 0.1);
 sample_num = xrow;
caseNum = sum(state);
controlNum = sample_num - caseNum;
threshold = caseNum/controlNum;

Hcase = 0; Hcontrol = 0;
Lcase = 0; Lcontrol = 0;


for i = 1:length(Data(:,1))
    if ProCaseToControl(i) > threshold
        Hcase = Hcase + disease(i);
        Hcontrol = Hcontrol + control(i);
    else
        Lcase = Lcase + disease(i);
        Lcontrol = Lcontrol + control(i);
    end
end
   
TP = Hcase; FP = Hcontrol;  H = TP + FP;
FN = Lcase; TN = Lcontrol;   L = TN + FN;
T  = TP + TN;  F = FP + FN; ALL = H + L;
% 列联表，O观察值，E期望值
O = [TP, FP, H; 
    TN, FN, L; 
    T, F, ALL];
E = [H*T/ALL, F*H/ALL;
      L*T/ALL, F*L/ALL]; 
err_rate = 0.5*( FN/(TP + FN) + FP / (FP + TN)); 
LR = 0;
for i=1:2
    for j=1:2
        LR = LR + O(i,j)*log(O(i,j)/E(i,j));
    end
end
LR = 2*LR;

% PSI (Predictive Summary Index) (Linn and Grunau 2006)
PSI = TP / (TP + FP) + TN / (TN + FN) - 1;

% Classification Error (CE)
  CE = (FP + FN)/(FP + FN + TP + TN);
 % Edclidean Distance (ED)
   ED = sqrt( (TP/(TP+FN) - 1)^2 + (TP/(TP + FP) - 1)^2 );
   
%  F-Measure
   beta = 0.5;
    F = (beta^2 + 1) * (TP/(TP + FP) * (TP/(TP + FN))) / ( beta^2 * TN/(TN+FN) * TN/(TN+FP));
 
  CA = 1 - CE;


end