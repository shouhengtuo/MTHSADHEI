%Tuo SH et al. MTHSA-DHEI:MTHSA-DHEI: Multitasking Harmony Search Algorithm for Detecting High-Order Epistatic Interactions
 clear;


warning('off');

%% harmony search algorithm parameters setting
for  F = 20:2:20
         PAR=0.5;
      
         HMCR = 0.98;
         TP0 = 0.5;
       Dim =   100;
       sample_num = 3000;
    %%
    folder = ['ResultData_PAR0.5TP0.5F=',num2str(F)];
    if ~exist(folder)
       mkdir(folder);
    end
    folder =[folder,'\'];

    %     A = {'Dim','Power','meanFEs','MeanTime', 'succFEs', 'succTime'};
    %     sheet = 1;
    %    xlRange = 'b1';
    %    xlswrite(dataFile,A,sheet,xlRange)

    for dataIndex =[10:12]
            switch(dataIndex)
                 case 1
                    epi_dim = 3;
                    startIndex = 1; endIndex = 100;
                    filepath = '.\threewayBests\'; filename = 'threewayBests'; 
                case 2
                    epi_dim = 3;
                    startIndex = 401; endIndex = 500;
                    filepath = '.\HWthreewayBests\';filename = 'HWthreewayBests'; 
                case 3
                    epi_dim = 4;
                    startIndex = 101; endIndex = 200; 
                    filepath = '.\fourwayBests\';filename = 'fourwayBests';
                case 4
                    epi_dim = 4;
                    startIndex = 1101; endIndex = 1200;
                    filepath = '.\fourwayNoLowBests\';    filename = 'fourwayNoLowBests';


%               case 9 %additive models
%                     epi_dim = 5;  Dim = 500;
%                     startIndex =001; endIndex = 10;
%                     filepath = '.\additive_model\5_0.10_0.10_EDM-1\';
%                     filename = '5_0.10_0.10_EDM-1';
% 
%              case 10 %additive models
%                     epi_dim = 5; Dim = 500;
%                     startIndex =001; endIndex = 10;
%                     filepath = '.\additive_model\5_0.25_0.10_EDM-1\';
%                     filename = '5_0.25_0.10_EDM-1';
% 
% 
%              case 13 %threshold models2
%                     epi_dim = 5;  Dim = 500;
%                     startIndex =001; endIndex = 10;
%                     filepath = '.\threshold_model\5_0.10_0.10_EDM-1\';
%                     filename = '5_0.10_0.10_EDM-1';
% 
%              case 14 %threshold models2
%                     epi_dim = 5;
%                     startIndex =001; endIndex = 10;
%                     filepath = '.\threshold_model\5_0.10_0.25_EDM-1\';
%                     filename = '5_0.10_0.25_EDM-1';
                      case 10 %additive models
                            epi_dim = 5; Dim = 500;
                            startIndex =001; endIndex = 100;
                            filepath = 'F:\diseaseModels\additive_model\5_0.25_0.10_EDM-1\';
                            filename = '5_0.25_0.10_EDM-1';

                     case 11 %threshold models2
                            epi_dim = 5;  Dim = 500;
                            startIndex =001; endIndex = 100;
                            filepath = 'F:\diseaseModels\threshold_model\5_0.10_0.10_EDM-1\';
                            filename = '5_0.10_0.10_EDM-1';


                      case 12 %multiplicative models
                            epi_dim = 4;  Dim = 1000;
                            startIndex =001; endIndex = 100;
                            filepath = 'F:\diseaseModels\Dim4_SNP1000_Samples4000\4_0.1_0.1_EDM-1\';
                            filename = '4_0.1_0.1_EDM-1'; 
                
        


            end

            dataNum = endIndex - startIndex +1;

           %% harmony search algorithm parameters setting  
               clsThreshold = [0.95,0.85, 0.75, 0.7,0.65,0.6,0.55, 0.53, 0.51, 0.5];
            M = length(clsThreshold);
            Ac = zeros(1,M);
             TP = zeros(1,M);
            FP = zeros(1,M);
            TN = zeros(1,M);
            FN = zeros(1,M);
            power = 0;
            power2 = 0;
            power3 = zeros(1,M);
            mdrPower = zeros(1,M);
            mdrscore = 0;

            %treeThreshold = 0.7;
            Evalutation_Times = [];
            Dim_FEs = [];
            TIME = [];
             succ = 0; 
            succEvalutation_Times = [];
            succTIME = [];


            for dataSet = startIndex:endIndex
               % a = dlmread(strcat(filepath,strcat('best',num2str(dataSet),'.txt')),',',1,0);
              if dataIndex>8 && dataIndex <=11
                  if dataSet < 10
                      ss = strcat('00',num2str(dataSet));
                  elseif dataSet < 100
                      ss = strcat('0',num2str(dataSet));
                  else
                      ss = '100';
                  end
                 data = dlmread([filepath,filename,'_',ss,'.txt'],'\t',1,0);

              elseif dataIndex>=12 && dataIndex <=20
                  if dataSet < 10
                      ss = strcat('00',num2str(dataSet));
                  elseif dataSet < 100
                      ss = strcat('0',num2str(dataSet));
                  else
                      ss = '100';
                  end
                 data = dlmread([filepath,ss,'.txt'],'\t',1,0);
              else
               a = dlmread(strcat(filepath,strcat('best',num2str(dataSet),'.txt')),'\t',1,0);
               AA = 0; Aa = 0; aa = 0;
                    for i =1:1500 %1:sample_num
                        for j = 1:epi_dim
                            if a(i,j) == 2
                                aa = aa + 1;
                            elseif a(i,j) == 1
                                Aa = Aa + 1;
                            else
                                AA = AA + 1;
                            end
                        end
                    end
                    AA = 2*AA /(sample_num*epi_dim);
                    Aa = 2*Aa /(sample_num*epi_dim);
                    aa = 2*aa /(sample_num*epi_dim);

                    b = zeros(sample_num,Dim-epi_dim);
                    for i = 1 : sample_num
                        for j = 1:Dim-epi_dim
                            r = rand;
                            if r <= aa
                                b(i,j) = 2;
                            elseif r <= aa+Aa
                                b(i,j) = 1;
                            else
                                b(i,j) = 0;
                            end
                              % b(i,j) = fix(rand*3);
                        end
                    end
                    data = [b,a];
              end

              %% parameter settings
              [m,n] = size(data);
              Dim = n - 1;



           max_iter = min(50*epi_dim*Dim, 4*epi_dim*50000) ;
             HMS =2*max(100, epi_dim*min(Dim/10,100)); %epi_dim*20;

            CX = Dim-epi_dim+1:Dim;




            p_value0 = max(1e-7,0.05/nchoosek(Dim,epi_dim));


        %% Search k-snp loci using Harmony search algorithm        
               tic;

                s = 2;
               [Task,NC,flag,Epi_Dim_FEs] = HS_2021_multiTask_UnifiedCoding2021(data, epi_dim, s, HMS,max_iter,CX, TP0, PAR,F,HMCR);

    %          [Candidate,canSize,Nc,runtime,flag] = HS_FOR_multiLOCI12(data,epi_dim,HMS,max_iter,maxIterForLocalSearch,CandidateSize,CX);
                 Evalutation_Times = [Evalutation_Times, NC];
                 Dim_FEs = [Dim_FEs,Epi_Dim_FEs];


                     runtime = toc;
                 TIME = [TIME  runtime];

               %% 1st stage power  
                if flag > 0
                    power = power + 1;
                     succ = succ + 1;
                    succEvalutation_Times = [succEvalutation_Times, NC];
                     succTIME = [succTIME  runtime];
                     fprintf('Index:%d,dataSet %3d: success search time(%f),   successTimes  %2d/%2d, FEs = %d\n',dataIndex,dataSet, runtime,succ, dataSet-startIndex+1,Epi_Dim_FEs);
                else 
                    fprintf(' Index:%d,dataSet %3d:  search time(%f) **** fail! *** \n',dataIndex,dataSet,runtime);
                end
             %% 2nd and 3rd stage 
                if flag > 0 && Gtest_score(data(:,CX),data(:,end)) < p_value0
                    power2 = power2 + 1;

                    [CE, CA, PSI,F_score] = MDR_2020_2(data(:,CX),data(:,end));
                  % TreeAcc = fitnessTreeBagger(data(:,CX),data(:,end));
                    mdrscore = mdrscore + CA;

                    for m = 1:length(clsThreshold)
                        if CA >clsThreshold(m)
                            mdrPower(1,m) = mdrPower(1,m) + 1;
                        end
    %                     if TreeAcc > clsThreshold(m)
    %                         treePower(1,m)  =treePower(1,m)  +  1;
    %                     end
                        if CA > clsThreshold(m) %&& TreeAcc > clsThreshold(m)
                            power3(1,m)  = power3(1,m) + 1;
                            TP(1,m) = TP(1,m)  + 1;

                        else
                            FN(1,m)  = FN(1,m)  + 1;
                        end
                        TN(1,m) = TN(1,m)  + 1;
                    end
                elseif flag == 0
                    FN(1,:)  = FN(1,:)  + 1;
                            for k = 1:4
                               X = Task(epi_dim-1).Criteria(k).bestXs;
                               L = length(X(:,1));
                                       for j = 1:L
                                           x = X(j,:);
                                                   if Gtest_score(data(:,x),data(:,end)) > p_value0
                                                       TN(1,:) = TN(1,:) + 1;
                                                       [CE, CA, PSI,F_score] = MDR_2020_2(data(:,x),data(:,end));
                                                       % TreeAcc = fitnessTreeBagger(data(:,x),data(:,end));

                                                             for m = 1:length(clsThreshold)
                                                                    if CA < clsThreshold(m) %|| TreeAcc < clsThreshold(m)
                                                                        TN(1,m) = TN(1,m) + 1;
                                                                    else
                                                                        FP(1,m) = FP(1,m) + 1;
                                                                    end
                                                             end
                                                   end
                                       end
                            end
                end

            end
         FN = dataNum - TP;
           meanFEs = mean(Evalutation_Times);
           meanTime = mean(TIME);
           Dim_FEs = mean(Dim_FEs);
           succFEs = mean(succEvalutation_Times);
           succTime = mean(succTIME);
           precision = TP./ (TP + FP);
           recall = TP./ (TP + FN);
           F_score = (2*precision .* recall)./(precision + recall);
           FDR = FP./(FP + TP);
           sensibility = TP./(TP + FN);
           specificity = TN./(TN + FP);

    %        % plot Roc curve
    %        figure(dataIndex);
    %        plot(1-specificity, sensibility,'*-');
    %        xlabel('False Positive Rate (1-Specificity)');
    %        ylabel('True Positive Rate(Sensitivity)');
    %        if dataIndex <= 8
    %           title(['EINME',num2str(dataIndex),': ROC curve']);
    %        else
    %            title(['EIME',num2str(dataIndex-8),': ROC curve']);
    %        end

    %        savefig(['.\roc_figures\', num2str(dataIndex),'.fig']);


            Results = [Dim,power, meanFEs, meanTime, succFEs, Dim_FEs,succTime, power, power2, power3, mdrPower,TP,FP,TN,FN, precision,recall,F_score,sensibility,specificity, FDR,mdrscore];
    %        dataFile = [folder,'Results.xls'];
          dataFile = strcat(folder,'MultiPower_',num2str(Dim));
             sheet = 1;
           xlRange = strcat('B', num2str(dataIndex+1)) ;
           xlswrite(dataFile,Results,sheet,xlRange)
            sheet2 = 2;
           range=strcat(['A'+2*(dataIndex-1),'2']);
           xlswrite(dataFile,[Evalutation_Times',TIME'],sheet2,range);
    end
end



