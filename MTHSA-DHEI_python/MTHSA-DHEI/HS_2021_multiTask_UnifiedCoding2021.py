import random
import math

from multi_criteriaEvaluationFuns2021 import *


class Task(object):
    pass


class Elite:
    pass


# Initialize multitasking Harmony Memories
def initial_Task(HMS, epi_dim, K, SNPs, EliteSize, data, state, FitNum):
    task = []
    X = np.zeros([K, HMS, epi_dim], dtype=int)
    HM = np.zeros([K, HMS, epi_dim], dtype=int)
    Fit = np.zeros([K, HMS, FitNum])
    eliteX = np.zeros([K, FitNum, EliteSize, epi_dim], dtype=int)
    eliteHM = np.zeros([K, FitNum, EliteSize, epi_dim], dtype=int)
    eliteFit = np.zeros([K, FitNum, EliteSize, FitNum])
    maxFit = np.zeros([K, FitNum])
    for i in range(K):
        task.append(Task())
        for j in range(HMS):
            SNP = random.sample(range(1, SNPs), epi_dim)
            SNP2 = np.sort(SNP[:i + 2])
            X[i, j, :] = np.concatenate((SNP2, SNP[i + 2:]), axis=0)
            HM[i, j, :] = SNP
            Fit[i, j, :] = multi_criteriaEvaluationFuns(data[:, SNP2], state)
        setattr(task[i], 'X', X[i, :, :])  # Store ordered solutions in X
        setattr(task[i], 'HM', HM[i, :, :])  # Store unordered solutions in HM
        setattr(task[i], 'Fit', Fit[i, :, :]) # the association score between SNP combination and disease status
        # Store the elite solutions in Elite
        elite = []
        for f in range(FitNum):
            elite.append(Elite())
            eliteX[i, f, :, :] = X[i, :EliteSize, :]
            eliteHM[i, f, :, :] = HM[i, :EliteSize, :]
            eliteFit[i, f, :, :] = Fit[i, :EliteSize, :]
            for e in range(EliteSize, HMS):
                max_Fit_inex = np.argmax(eliteX[i, f, :, f])
                if Fit[i, e, f] < eliteX[i, f, max_Fit_inex, f]:
                    eliteX[i, f, max_Fit_inex, :] = X[i, e, :]
                    eliteHM[i, f, max_Fit_inex, :] = HM[i, e, :]
                    eliteFit[i, f, max_Fit_inex, :] = Fit[i, e, :]
            setattr(elite[f], 'X', eliteX[i, f, :, :])
            setattr(elite[f], 'HM', eliteHM[i, f, :, :])
            setattr(elite[f], 'Fit', eliteFit[i, f, :, :])
        setattr(task[i], 'Elite', elite)
        maxFit[i, :] = np.max(Fit[i, :, :], axis=0)
    maxFit = np.max(maxFit, axis=0)
    return task, maxFit


#  Improvise a new SNP in Task_k
def GenerateNewHamronyInTask_k(task, EliteSize, epi_dim, HMCR, PAR , HMS, k, d, d0, bestNum, i, F, SNPs):
    if random.uniform(0, 1) < HMCR:
        a = math.ceil(random.uniform(0, 1) * HMS) - 1
        b = math.ceil(random.uniform(0, 1) * epi_dim) - 1
        Xnewi = task[k].X[a, b]
        if random.uniform(0, 1) < PAR:
            sPar = math.ceil(random.uniform(0, 1) * 4)
            b = math.ceil(random.uniform(0, 1) * epi_dim) - 1
            c = math.ceil(random.uniform(0, 1) * EliteSize) - 1
            # while c == a:
            #     c = math.ceil(random.uniform(0, 1) * EliteSize) - 1
            bs = math.ceil(random.uniform(0, 1) * bestNum) - 1
            # three kinds of strategies
            if sPar == 1:
                Xnewi = task[k].Elite[d].Xbest[bs, b]
            elif sPar == 2:
                Xnewi = task[k].Elite[d].X[bs, b]
            elif sPar == 3:
                e = math.ceil(random.uniform(0, 1) * HMS) - 1
                L = task[k].Elite[d].X[c, b] - task[k].X[e, i]
                Xnewi = round(Xnewi + F * random.uniform(0, 1) * L)
                # Boundary detection
                if Xnewi > SNPs:
                    Xnewi = SNPs - max(0, round(random.normalvariate(0, min(10, max(L / 10, 1)))))
                elif Xnewi < 1:
                    Xnewi = 1 + max(0, round(random.normalvariate(0, min(10, max(L / 10, 1)))))
            else:
                e = math.ceil(random.uniform(0, 1) * HMS) - 1
                L = task[k].Elite[d].Xbest[bs, b] - task[k].X[e, b]
                Xnewi = round(Xnewi + F * random.uniform(0, 1) * L)
                if Xnewi > SNPs:
                    Xnewi = SNPs - max(0, round(random.normalvariate(0, min(10, max(L / 10, 1)))))
                elif Xnewi < 1:
                    Xnewi = 1 + max(0, round(random.normalvariate(0, min(10, max(L / 10, 1)))))
    else:
        Xnewi = np.ceil(random.uniform(0, 1) * SNPs) - 1
    return Xnewi

#  Improvise a new SNP from Task_k for knowledge transfer
def GenerateNewHamronyFromTask_k(task, EliteSize, epi_dim, HMCR, PAR , HMS, k, d, d0, bestNum, i, F, SNPs):
    if random.uniform(0, 1) < HMCR:
        a = math.ceil(random.uniform(0, 1) * EliteSize) - 1
        b = math.ceil(random.uniform(0, 1) * epi_dim) - 1
        Xnewi = task[k].Elite[d].X[a, b]  # select SNP from d-th Elite 
        if random.uniform(0, 1) < PAR:
            sPar = math.ceil(random.uniform(0, 1) * 3)
            b = math.ceil(random.uniform(0, 1) * epi_dim) - 1
            c = math.ceil(random.uniform(0, 1) * EliteSize) - 1
            while c == a:
                c = math.ceil(random.uniform(0, 1) * EliteSize) - 1
            bs = math.ceil(random.uniform(0, 1) * bestNum) - 1
            # three kinds of strategies
            if sPar == 1:
                Xnewi = task[k].Elite[d].Xbest[bs, b]
            elif sPar == 2:
                Xnewi = task[k].Elite[d].X[bs, b]
           
            else:
                e = math.ceil(random.uniform(0, 1) * HMS) - 1
               
                
                L = task[k].Elite[d0].Xbest[bs, b] - task[k].X[e, b]
                Xnewi = round(Xnewi + F * random.uniform(0, 1) * L)
                if Xnewi > SNPs:
                    Xnewi = SNPs - max(0, round(random.normalvariate(0, min(10, max(L / 10, 1)))))
                elif Xnewi < 1:
                    Xnewi = 1 + max(0, round(random.normalvariate(0, min(10, max(L / 10, 1)))))
    else:
        Xnewi = np.ceil(random.uniform(0, 1) * SNPs) - 1
    return Xnewi

def Array_BcontiansA(ArrayA, ArrayB):
    ArrayA = np.sort(ArrayA)
    ArrayB = np.sort(ArrayB)
    a = np.size(np.unique(ArrayB, axis=0), axis=0)
    b = np.size(np.unique(np.append([ArrayA], ArrayB, axis=0), axis=0), axis=0)
    if a == b:
        return True
    else:
        return False


#  Adjust Xnew
def XnewAdjust(Xnew, task, dim, SNPs, k, FitNum ):
    Xtemp = Xnew
    Xnew0 = np.sort(Xnew[:dim])
    Xnew[:dim] = Xnew0
    flag_Xnew_dim = Array_BcontiansA(Xnew[:dim], task[k].X[:, :dim])
    if flag_Xnew_dim == True:
        j = math.ceil(random.uniform(0, 1) * dim) - 1
        r = math.ceil(random.uniform(0, 1) * SNPs) - 1
        while r in Xnew:
            r = math.ceil(random.uniform(0, 1) * SNPs) - 1
        Xnew[j] = r
        Xtemp = Xnew
        Xnew0 = np.sort(Xnew[:dim])
        Xnew[:dim] = Xnew0
    for b in range(FitNum):
        Xtemp = Xnew
        Xnew0 = np.sort(Xnew[:dim])
        Xnew[:dim] = Xnew0
        flag_Xnew_dim = Array_BcontiansA(Xnew[:dim], task[k].Elite[b].X[:, :dim])
        if flag_Xnew_dim == True:
            j = math.ceil(random.uniform(0, 1) * dim) - 1
            r = math.ceil(random.uniform(0, 1) * SNPs) - 1
            while r in Xnew:
                r = math.ceil(random.uniform(0, 1) * SNPs) - 1
            Xnew[j] = r
            Xtemp = Xnew
            Xnew0 = np.sort(Xnew[:dim])
            Xnew[:dim] = Xnew0
    return Xnew, Xtemp


# Improvise a new harmony Xnew
def GenerateNewHarmony(task, SNPs, EliteSize, epi_dim, HMCR, PAR, HMS, k, FitNum, dim, Ds, Rs, TP, bestNum, K, F):
    Xnew = np.zeros(epi_dim, dtype=int)
    Xtemp = np.zeros(epi_dim, dtype=int)
    d = math.ceil(random.uniform(0, 1) * FitNum) - 1
    i = 0
    while i <= epi_dim-1:
        k0 = math.ceil(random.uniform(0, 1) * K) - 1
        while k0 == k:
            k0 = math.ceil(random.uniform(0, 1) * K) - 1
        if Rs >= TP:
            # Generate a new harmony from the memories of the current task
            d0 = d
            Xnew[i] = GenerateNewHamronyInTask_k(task, EliteSize, epi_dim, HMCR, PAR, HMS, k, d, d0, bestNum, i, F, SNPs)
        else:
            # Generate a new harmony by transferring learning from another task
            d0 = math.ceil(random.uniform(0, 1) * FitNum) - 1
            while d0 == d:
                d0 = math.ceil(random.uniform(0, 1) * FitNum) - 1
            Xnew[i] = GenerateNewHamronyFromTask_k(task, EliteSize, epi_dim, HMCR, PAR, HMS, k0, d, d0, bestNum, i, F, SNPs)
        if i == 0 or Xnew[i] not in Xnew[:i]:
            i = i + 1
            Xnew, Xtemp = XnewAdjust(Xnew, task, dim, SNPs, k, FitNum)
        if i == epi_dim:
            if Rs < TP:
                Xnew, Xtemp = XnewAdjust(Xnew, task, Ds, SNPs, k0, FitNum)
    Xnew = Xnew.astype(int)
    return Xnew, Xtemp


# Updating the Harmony Memory and elite Harmony Memory of task k
def Update_HM(task, EliteSize, epi_dim, HMS, FitNum, dim, NC, Epi_Dim_FEs, Xnew, XnewScore, k, maxFit, Xtemp, Epi_Dim, bestNum, max_iter, data, state):
    for i in range(HMS):
        sn = []
        for j in XnewScore[:(FitNum - 1)]:
            sna = [x for x in task[k].Fit[i, :(FitNum - 1)] if x > j]
            sn = np.append(sna, sn)
        if np.size(sn) > 2 or XnewScore[FitNum - 1] < task[k].Fit[i, (FitNum - 1)] and random.uniform(0, 1) < (
                1 - NC / max_iter):
            task[k].X[i, :] = Xnew
            task[k].HM[i, :] = Xtemp
            task[k].Fit[i, :] = XnewScore
        break
    for i in range(FitNum):
        fworst = max(task[k].Elite[i].Fit[:, i])
        worstId = np.argmax(task[k].Elite[i].Fit[:, i])
        if fworst > XnewScore[i]:
            task[k].Elite[i].X[worstId, :] = Xnew
            task[k].Elite[i].HM[worstId, :] = Xtemp
            task[k].Elite[i].Fit[worstId, :] = XnewScore
            for s in range(bestNum):
                if XnewScore[i] < task[k].Elite[i].fbest[s]:
                    task[k].Elite[i].fbest[s] = XnewScore[i]
                    task[k].Elite[i].Xbest[s] = Xnew
                    if dim < epi_dim:
                        Score = multi_criteriaEvaluationFuns(data[:, Xnew[:(dim + 2)]], state)
                        Score = Score / maxFit
                        if dim + 1 == Epi_Dim:
                            Epi_Dim_FEs = Epi_Dim_FEs + 1
                        NC = NC + 1
                        for si in range(FitNum):
                            for sj in range(EliteSize):
                                if Score[si] < task[dim-1].Elite[si].Fit[sj, si]:
                                    task[dim-1].Elite[si].X[sj, :] = np.sort(Xnew)
                                    task[dim-1].Elite[si].HM[sj, :] = Xnew
                                    task[dim-1].Elite[si].Fit[sj, :] = Score
                                    break
                    break
    return task, NC, Epi_Dim_FEs


#  the main function of MTHSA-DHEI
def MTHSA(sample_data, epi_dim, HMS, s, max_iter, CX):
    data = sample_data[:, :-1]
    state = sample_data[:, -1]
    # initial arguments
    Epi_Dim = epi_dim  # the dimension of functional snp combination
    epi_dim = epi_dim + s  # tThe highest order of SNP interactions to be detected
    HMCR = 0.98  # Harmony memory considering rate
    PAR = 0.5  # Pitching adjusting rate
    TP = 0.35  # Transferring learning rate
    F = 5
    fdim = len(CX)  # the dimension of functional SNP combination
    n = data.shape[1]
    SNPs = n - 1  # number of SNPs in sample dataset
    K = epi_dim - 1  # The number of tasks for detecting high-order SNP epistatic interactions
    FitNum = 4   # 4 evaluation criteria that MTSHA-DHEI employed
    bestNum = epi_dim
    EliteSize = min(5 * epi_dim, HMS)  # the number of elite solutions in elite
    # Parameter for recording the algorithm status
    NC = 0
    Epi_Dim_FEs = 0
    flag = -1
    # initial tasks
    task, maxFit = initial_Task(HMS, epi_dim, K, SNPs, EliteSize, data, state, FitNum)
    NC = NC + HMS * K
    Epi_Dim_FEs = Epi_Dim_FEs + HMS
    # Normalize all of fitness values
    for i in range(K):
        task[i].Fit = task[i].Fit / maxFit
        for j in range(FitNum):
            task[i].Elite[j].Fit = task[i].Elite[j].Fit / maxFit
            bestId = np.argsort(task[i].Elite[j].Fit[:, j])
            X = task[i].Elite[j].X[bestId, :]
            Fit = task[i].Elite[j].Fit[bestId, j]
            task[i].Elite[j].Xbest = X[:bestNum,:]
            task[i].Elite[j].fbest = Fit[:bestNum]
    # HM　search
    # max_iter used as the terminal condition for HS.
    while NC < max_iter:
        for dim in range(2, epi_dim+1):  # the dimensionality of each task
            k = dim - 2  # index of tasks in HM
            Rs = random.uniform(0, 1)
            Ks = k
            Ds = dim
            # Improvise a new harmony
            Xnew, Xtemp  = GenerateNewHarmony(task, SNPs, EliteSize, epi_dim, HMCR, PAR, HMS, k, FitNum, dim, Ds, Rs, TP, bestNum, K, F)
            # Update HM and Elite-HM of each task
            if Rs < TP:
                XnewScore = multi_criteriaEvaluationFuns(data[:, Xnew[:(Ds + 1)]], state)
                XnewScore = XnewScore / maxFit
                if Ds == Epi_Dim:
                    Epi_Dim_FEs = Epi_Dim_FEs + 1
                NC = NC + 1
                task, NC, Epi_Dim_FEs = Update_HM(task, EliteSize, epi_dim, HMS, FitNum, dim, NC, Epi_Dim_FEs, Xnew, XnewScore, Ks, maxFit,
                          Xtemp, Epi_Dim, bestNum, max_iter, data, state)
            else:
                XnewScore = multi_criteriaEvaluationFuns(data[:, Xnew[:(dim + 1)]], state)
                XnewScore = XnewScore / maxFit
                if dim == Epi_Dim:
                    Epi_Dim_FEs = Epi_Dim_FEs + 1
                NC = NC + 1
                task, NC, Epi_Dim_FEs = Update_HM(task, EliteSize, epi_dim, HMS, FitNum, dim, NC, Epi_Dim_FEs, Xnew, XnewScore, k, maxFit,
                          Xtemp, Epi_Dim, bestNum, max_iter, data, state)
            # Check the termination Conditions
            cflag = 0
            for ci in range(fdim):
                if CX[ci] in Xnew:
                    cflag = cflag + 1
            if cflag == fdim:
                task[fdim - 1].Elite[1].X[1, :] = Xnew
                task[fdim - 1].Elite[1].Fit[1, :] = XnewScore
                flag = 1
                break
        if flag == 1:
            break
    return Task, NC, flag, Epi_Dim_FEs


# filepath = "./NDME-1/"
# data = np.loadtxt(filepath +'001.txt', delimiter='\t')
# epi_dim = 3
# s = 2
# HMS = 50
# max_iter = 20000
# CX = [97,98,99]
# task, NC, flag, Epi_Dim_FEs = MTHSA(data, epi_dim, HMS, s, max_iter, CX)
# print(NC, flag, Epi_Dim_FEs)