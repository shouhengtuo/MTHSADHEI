 function [Task,NC,flag,Epi_Dim_FEs] = HS_2021_multiTask_UnifiedCoding2021(data,epi_dim,s,HMS,max_iter,CX, TP, PAR,F,HMCR)
% 多任务统一编码
%input--------------------------------------------------------------------
% data-----------------input dataset
% epi_dim--------------the epistasis dimension
% HMS--------------the size of harmony memory(HM)

%%-------------------------------------------------------------------------
% initial arguments
Epi_Dim = epi_dim;
Epi_Dim_FEs = HMS;
%F = 5;
epi_dim = epi_dim + s;
% HMCR=0.98;
% PAR=0.7;
% TP = 0.35;
fdim = length(CX);
n=size(data,2);
State=data(:,n);
K = epi_dim-1 ; % 任务数 2-->epi_dim 
FitNum = 4;
bestNum = epi_dim;
%% ---------------------------------------------------------------
EliteSize =min(10*epi_dim,ceil(HMS/5));
% CX
NC=0;
SNPs=n-1;  %% 总SNP个数
flag = -1;
maxFit = [];
for k = 1:K
    %% 初始化
    
    dim = k+1;
    Task(k).X=zeros(HMS,epi_dim);
    snp=[];
    for i=1:HMS
        snp(1)=ceil(rand*SNPs);
        for j=2:epi_dim
          snp(j)=ceil(rand*SNPs); 
          while ismember(snp(j),snp(1:j-1)) 
             snp(j)=ceil(rand*SNPs);        
          end
        end
        temp=snp;
        snp2=sort(snp(1:(k+1)));
        while ismember(snp,Task(k).X,'rows')
            j=ceil(rand*epi_dim);
            snp(j)=ceil(rand*SNPs); 
            temp=snp;
            snp2=sort(snp(1:(k+1)));
        end

        X(i,:)=[snp2,snp(k+2 : epi_dim)];   %% X中存放有序的解
        HM(i,:)=temp;  %% HM中相应存放无序解
       % K2Score,GtestP_value,Gini_Score,JE_Score
        [Fit(i,1),Fit(i,2),Fit(i,3),Fit(i,4)] = multi_criteriaEvaluationFuns2021(data(:,HM(i,1:dim)),State);
        NC = NC + 1;
        snp=[];
    end
    maxFit = [maxFit; Fit];
       
    Task(k).X = X(1:HMS,:);
    Task(k).HM = HM(1:HMS,:);
    Task(k).Fit = Fit(1:HMS,:);
    
    for i = 1:FitNum
        Elite(i).X = X(1:EliteSize,:);
        Elite(i).HM = HM(1:EliteSize,:);
        Elite(i).Fit = Fit(1:EliteSize,:);
    end
    
        for j = 1:EliteSize
            for s = 1:FitNum
                for i = EliteSize+1 : length(Fit(:,1))
                    if Elite(s).Fit(j,s) > Fit(i,s)
                        Elite(s).Fit(j,:) = Fit(i,:);
                        Elite(s).X(j,:) = X(i,:);
                        Elite(s).HM(j,:) = HM(i,:);
                        break;
                    end
                end
            end          
        end
    
    Task(k).Elite = Elite;
   
   
end

maxFit = max(maxFit);
% 对所有评价指标适应值进行归一化
for k = 1:K
    Task(k).Fit = Task(k).Fit ./ maxFit;
    for i = 1:FitNum
        Task(k).Elite(i).Fit = Task(k).Elite(i).Fit ./maxFit;      
        
        [Task(k).Elite(i).fbest, bestId] = mink(Task(k).Elite(i).Fit(:,i),bestNum);
        Task(k).Elite(i).Xbest = Task(k).Elite(i).X(bestId,:);
    end
end


[s(1),s(2),s(3),s(4)] =  multi_criteriaEvaluationFuns2021(data(:,CX),State); 
s = s./maxFit;
 
Dims = [2:epi_dim]; 

 LT=0;
%%-------------------------------------------------------------------------
tic;
while NC <= max_iter 
   for dim = Dims
       k = dim - 1;  %% 从第k个任务群众探索     
      Rs = rand;
      
      if Rs < TP %% 迁移学习
          k0 = ceil(rand*K);
          while k0 == k
             k0 = ceil(rand*K);
          end
      end
   
        
             Ks = k;
             Ds = dim;
       
      i=1;
       d = ceil(rand*FitNum);
      while i<= epi_dim
         
            if Rs >= TP %% 从当前任务中进行 优化组合
                 if rand<HMCR
                     a = ceil(rand*EliteSize);                   
                     b = ceil(rand*epi_dim);
                     Xnew(i) = Task(k).X(a,b);
                     if rand < PAR
                         sPar = ceil(rand*4);
                         b = ceil(rand*epi_dim);
                         c = ceil(rand*EliteSize); 
                         while c == a
                             c = ceil(rand*EliteSize); 
                         end
                          bs = ceil(rand*bestNum);
                         switch sPar
                             case 1
                                 Xnew(i) = Task(k).Elite(d).Xbest(bs,b);
                             case 2
                                  e = ceil(rand*HMS);
                                  
                                  L = Task(k).Elite(d).X(c,b)-Task(k).X(e,i);
                                  
                                 % Xnew(i) = round(Xnew(i) + 0.5 * (Task(k).Elite(d).X(c,b)-Task(k).X(e,i)));
                                 Xnew(i) = round( Xnew(i) + F * rand * L);
                                      
                                 if Xnew(i) > SNPs || Xnew(i) <1 
                                     Xnew(i) = ceil(rand*SNPs);
                                 end
                             case 3   
                                 e = ceil(rand*HMS);
                                 L = Task(k).Elite(d).Xbest(bs,b) - Task(k).X(e,b);
                                 Xnew(i)=round(Xnew(i) + F * rand * L);
                                 %Xnew(i)=max(min(Xnew(i),SNPs),1);
                                         if  Xnew(i) > SNPs
                                             Xnew(i) = SNPs - max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                         elseif Xnew(i) < 1
                                             Xnew(i) = 1 + max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                         end
                         end
                     end
                     
                 else
                     Xnew(i)=ceil(rand*SNPs);
                 end
            else    
                 if rand<HMCR                     
                     a = ceil(rand*EliteSize);
                     b = ceil(rand*epi_dim);
                   
                     Xnew(i) = Task(k0).Elite(d).X(a,b);
                    if rand < PAR
                         sPar = ceil(rand*3);                      
                         b = ceil(rand*epi_dim);    
                         c = ceil(rand*EliteSize);
                         while c == a
                             c = ceil(rand*EliteSize);
                         end
                         bs = ceil(rand*bestNum);
                         switch sPar
                             case 1
                                 Xnew(i) = Task(k0).Elite(d).Xbest(bs,b);                            
                             case 2
                                 e = ceil(rand*HMS);
                                 L = Task(k0).Elite(d).X(bs,b)-Task(k0).X(e,i);
                                 Xnew(i) = Xnew(i)+round(rand * (Task(k0).Elite(d).X(bs,b)-Task(k0).X(e,i)));
                                 if Xnew(i) > SNPs || Xnew(i) <1 
                                     Xnew(i) = ceil(rand*SNPs);
                                 end
%                                      Xnew(i) = round( Xnew(i) + 2 * rand * L);
%                                         if  Xnew(i) > SNPs
%                                              Xnew(i) = SNPs - max(0,round(normrnd(0,min(10,max(L/10,1)))));
%                                          elseif Xnew(i) < 1
%                                              Xnew(i) = 1 + max(0,round(normrnd(0,min(10,max(L/10,1)))));
%                                          end
                                 

                             case 3
                                 d0 = ceil(rand*FitNum);
                                 while d0 == d
                                     d0 = ceil(rand*FitNum);
                                 end
                                 e = ceil(rand*HMS);
                                 L = Task(k0).Elite(d0).Xbest(bs,b) - Task(k0).X(e,b);
                                 Xnew(i) = Xnew(i)+round( F * rand * L);
                                 if Xnew(i) > SNPs || Xnew(i) <1 
                                     Xnew(i) = ceil(rand*SNPs);
                                 end
                                       if  Xnew(i) > SNPs
                                             Xnew(i) = SNPs - max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                         elseif Xnew(i) < 1
                                             Xnew(i) = 1 + max(0,round(normrnd(0,min(10,max(L/10,1)))));
                                         end
                         end
                     end
                 else
                        Xnew(i)=ceil(rand*SNPs);
                 end    
                 
            end
          if i==1 || ~ismember(Xnew(i),Xnew(1:i-1))
              i = i + 1;         
          end
                
          if ( i-1 == epi_dim )              
              Xtemp=Xnew;
              Xnew0=sort(Xnew(1:dim));
              Xnew(1:dim) = Xnew0;
                  while ( ismember(Xnew(1:dim),Task(k).X(:,1:dim),'rows')  )
                     j=ceil(rand*dim);
                      r=ceil(rand*SNPs);
                      while ismember(r,Xnew)
                          r=ceil(rand*SNPs);
                      end
                      Xnew(j)=r;
                      Xtemp=Xnew;
                      Xnew0=sort(Xnew(1:dim));                  
                      Xnew(1:dim) = Xnew0;
                  end
                  
                      for b = 1:FitNum
                          Xtemp=Xnew;
                          Xnew0=sort(Xnew(1:dim));
                          Xnew(1:dim) = Xnew0;
                          while ( ismember(Xnew(1:dim),Task(k).Elite(b).X(:,1:dim),'rows')  )
                              j=ceil(rand*dim);
                              r=ceil(rand*SNPs);
                              while ismember(r,Xnew)
                                  r=ceil(rand*SNPs);
                              end
                              Xnew(j)=r;
                              Xtemp=Xnew;
                              Xnew0=sort(Xnew(1:dim));                  
                              Xnew(1:dim) = Xnew0;
                          end
                      end
                if Rs < TP
                    Xtemp=Xnew;
                    Xnew0=sort(Xnew(1:Ds));
                    Xnew(1:Ds) = Xnew0;
                    while ( ismember(Xnew(1:Ds),Task(k0).X(:,1:Ds),'rows')  )
                         j=ceil(rand*Ds);
                          r=ceil(rand*SNPs);
                          while ismember(r,Xnew)
                              r=ceil(rand*SNPs);
                          end
                          Xnew(j)=r;
                          Xtemp=Xnew;
                          Xnew0=sort(Xnew(1:Ds));                  
                          Xnew(1:Ds) = Xnew0;
                    end
                  
                      for b = 1:FitNum
                          Xtemp=Xnew;
                          Xnew0=sort(Xnew(1:Ds));
                          Xnew(1:Ds) = Xnew0;
                          while ( ismember(Xnew(1:Ds),Task(k0).Elite(b).X(:,1:Ds),'rows')  )
                              j=ceil(rand*Ds);
                              r=ceil(rand*SNPs);
                              while ismember(r,Xnew)
                                  r=ceil(rand*SNPs);
                              end
                              Xnew(j)=r;
                              Xtemp=Xnew;
                              Xnew0=sort(Xnew(1:Ds));                  
                              Xnew(1:Ds) = Xnew0;
                          end
                      end
                end
          end  
      end
     
     
     
      if Rs < TP
             [XnewScore(1),XnewScore(2),XnewScore(3),XnewScore(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:Ds)),State); 
             XnewScore = XnewScore ./ maxFit;
             if Ds == Epi_Dim
                 Epi_Dim_FEs = Epi_Dim_FEs + 1;
             end
             NC=NC+1; 

             for i = 1:HMS        
                sn = find(XnewScore(1:(FitNum-1)) < Task(Ks).Fit(i,1:(FitNum-1)));
                if    length(sn) > 2 || (XnewScore(FitNum) < Task(Ks).Fit(i,FitNum) && rand < (1 - NC / max_iter))
                    %% 
                    Task(Ks).X(i,:) = Xnew;
                    Task(Ks).HM(i,:) = Xtemp;
                    Task(Ks).Fit(i,:) = XnewScore;
%                     fprintf('11\n');
                    break; % 只替换其中之一
                end
             end

            for i = 1:FitNum
                 [fworst,worstId] = max(Task(Ks).Elite(i).Fit(:,i));
                 if fworst > XnewScore(i)
                     Task(Ks).Elite(i).X(worstId,:) = Xnew;
                     Task(Ks).Elite(i).HM(worstId,:) = Xtemp;
                     Task(Ks).Elite(i).Fit(worstId,:) = XnewScore;
                         for s = 1:bestNum
                             if XnewScore(i) < Task(Ks).Elite(i).fbest(s)
                                 Task(Ks).Elite(i).fbest(s) = XnewScore(i);
                                 Task(Ks).Elite(i).Xbest(s,:) = Xnew;
                                  %% 扩展
                                    if dim < epi_dim
                                           [Score(1),Score(2),Score(3),Score(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:dim+1)),State); 
                                           Score = Score ./ maxFit;
                                            if dim+1 == Epi_Dim
                                                 Epi_Dim_FEs = Epi_Dim_FEs + 1;
                                             end
                                           NC = NC + 1;
                                           for si = 1:FitNum
                                               for sj = 1:EliteSize
                                                   if Task(dim).Elite(si).Fit(sj,si) > Score(si)
                                                       Task(dim).Elite(si).X(sj,:) = sort(Xnew);
                                                       Task(dim).Elite(si).HM(sj,:) = Xnew;
                                                       Task(dim).Elite(si).Fit(sj,:) = Score;
                                                       break;
                                                   end
                                               end
                                           end
                                    end
                                 
                                 
                                 break;
                             end
                         
                     end
                 end
            end
      else
          [XnewScore(1),XnewScore(2),XnewScore(3),XnewScore(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:dim)),State); 
         XnewScore  = XnewScore ./ maxFit;
             if dim == Epi_Dim
                 Epi_Dim_FEs = Epi_Dim_FEs + 1;
             end
              NC=NC+1; 
             for i = 1:HMS        
                sn = find(XnewScore(1:(FitNum-1)) < Task(k).Fit(i,1:(FitNum-1)));
                if    length(sn) >2 || (XnewScore(FitNum) < Task(k).Fit(i,FitNum) &&  rand < (1 - NC / max_iter))
                    %% 
                    Task(k).X(i,:) = Xnew;
                    Task(k).HM(i,:) = Xtemp;
                    Task(k).Fit(i,:) = XnewScore;
%                     fprintf('22\n');
                    break; % 只替换其中之一
                end
             end


                %% 更新 不同评价标准 的精英集合       
             for i = 1:FitNum
                [fworst,worstId] = max(Task(k).Elite(i).Fit(:,i));
                 if fworst > XnewScore(i)
                     Task(k).Elite(i).X(worstId,:) = Xnew;
                     Task(k).Elite(i).HM(worstId,:) = Xtemp;
                     Task(k).Elite(i).Fit(worstId,:) = XnewScore;
                         for s = 1:bestNum
                             if XnewScore(i) < Task(k).Elite(i).fbest(s)
                                 Task(k).Elite(i).fbest(s) = XnewScore(i);
                                 Task(k).Elite(i).Xbest(s,:) = Xnew;
                                 %% 扩展
                                    if dim < epi_dim
                                           [Score(1),Score(2),Score(3),Score(4)] = multi_criteriaEvaluationFuns2021(data(:,Xnew(1:dim+1)),State); 
                                           Score = Score ./ maxFit;
                                            if dim+1 == Epi_Dim
                                                 Epi_Dim_FEs = Epi_Dim_FEs + 1;
                                             end
                                           NC = NC + 1;
                                           for si = 1:FitNum
                                               for sj = 1:EliteSize
                                                   if Task(dim).Elite(si).Fit(sj,si) > Score(si)
                                                       Task(dim).Elite(si).X(sj,:) = sort(Xnew);
                                                       Task(dim).Elite(si).HM(sj,:) = Xnew;
                                                       Task(dim).Elite(si).Fit(sj,:) = Score;
                                                       break;
                                                   end
                                               end
                                           end
                                    end
                                 
                                 break;
                             end
                       
                     end                    
                 end
             end
      end
     %%  The program is terminted if the Xnew is the solution. 
     cflag = 0;
     for ci = 1:fdim
         
        if ismember(CX(ci), Xnew)
            
            cflag = cflag + 1;
        end
     end
     if cflag == fdim
         Task(fdim-1).Elite(1).X(1,:) = Xnew;
         Task(fdim-1).Elite(1).Fit(1,:) = XnewScore;
          flag = 1;            
          break;
      end
      

     
   end
  if flag == 1
      break;
  end
end

totaltime=toc;
% 