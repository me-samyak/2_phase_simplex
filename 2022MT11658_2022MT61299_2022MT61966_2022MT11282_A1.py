import numpy as np
def simplex_algo(filename="input.txt"):
    with open(filename, 'r') as file:
        sections = file.read().split('[')

    sections=sections[1:]
    for i in range(len(sections)):
        sections[i]='['+sections[i]

    A = np.array([list(map(float, line.split(','))) for line in sections[1][3:].strip().split('\n')])
    b = np.array(list(map(float, sections[2][3:].strip().split('\n'))))
    c = np.array(list(map(float, sections[4][3:].strip().split(','))))

    constraint_types = sections[3].strip().split('\n')
    constraint_types.pop(0)
    for i in range(len(constraint_types)-1):
        constraint_types[i]=constraint_types[i].lstrip()
        constraint_types[i]=constraint_types[i].rstrip()
    objective = sections[0].strip().split('\n')
    objective.pop(0)
    objective=objective[0]
    if (objective=="maximize"):
        c=c*(-1)
    for i in range(len(b)):
        if (b[i]<0):
            b[i]=b[i]*(-1)
            A[i]=A[i]*(-1)
            if (constraint_types[i]=='>='):
                constraint_types[i]='<='
            elif (constraint_types[i]=='<='):
                constraint_types[i]='>='
    count1=0
    count2=0
    count3=0
    for i in range(len(constraint_types)):
        if (constraint_types[i]=='='):
            count1+=1
        elif (constraint_types[i]=='>='):
            count2+=1
        elif (constraint_types[i]=='<='):
            count3+=1
    numzero=count1+2*(count2+count3)
    result = np.zeros((len(b),len(A[0])+numzero))
    modc=np.zeros(len(c)+count2+count3)
    for i in range(len(c)):
        modc[i]=c[i]
    for i in range(len(A)):
        for j in range(len(A[0])):
            result[i][j]=A[i][j]
    temp=len(A[0])
    for i in range(len(constraint_types)):
        if (constraint_types[i]=="="):
            continue
        elif (constraint_types[i]==">="):
            result[i][temp]=-1
            temp+=1
        else:
            result[i][temp]=1
            temp+=1
    for i in range(len(constraint_types)):
        result[i][temp]=1
        temp+=1
    initialtable=np.zeros((len(b),(len(result[0])+1)))
    for i in range(len(b)):
        initialtable[i][0]=b[i]
    for i in range(len(result)):
        for j in range(len(result[0])):
            initialtable[i][j+1]=result[i][j]
    numvar=len(result[0])
    numbasic=len(b)
    n=numvar-numbasic
    cost=(-1*np.sum(b))
    basicvarindex=np.zeros(numbasic)
    for i in range(numvar-numbasic+1,numvar+1):
        basicvarindex[i-(numvar-numbasic+1)]=i
    reduced_cost=np.zeros(numvar)
    for i in range(numvar - numbasic):
        reduced_cost[i] = -1 * np.sum(result[:, i]) 
    ############################################################################
    while True:
        result=np.round(result,6)
        reduced_cost=np.round(reduced_cost,6)
        cost=np.round(cost,6)
        b=np.round(b,6)
        if np.any(reduced_cost<0):
            pivotcol = np.where(reduced_cost < 0)[0][0]
            column_positive = np.any(result[:,pivotcol] > 0)
            if column_positive:
                mini=2**30-1
                for j in range(0,len(b)):
                    if result[:,pivotcol][j]>0:
                        mini=min(mini,b[j]/result[:,pivotcol][j])
                for j in range(0,len(b)):
                    if result[:,pivotcol][j]>0:
                        if b[j]/result[:,pivotcol][j]==mini:
                            index1=j
                            break
                for i in range(0,numbasic):
                    if (i!=index1):
                        alpha=-1*(result[:,pivotcol][i]/result[:,pivotcol][index1])
                        result[i]=result[i]+(alpha)*result[index1]
                        b[i]=b[i]+alpha*b[index1]
                alpha1=(-1*(reduced_cost[pivotcol])/result[:,pivotcol][index1])
                reduced_cost=reduced_cost+alpha1*result[index1]
                cost=cost+alpha1*b[index1]
                b[index1]=b[index1]/result[:,pivotcol][index1]
                result[index1]=result[index1]/result[:,pivotcol][index1]
                basicvarindex[index1]=pivotcol+1
            else:
                statusofsol="Infeasible"
                finaltable=np.zeros((len(b),len(result[0])+1))
                for i in range(len(b)):
                    finaltable[i][0]=b[i]
                for i in range(len(result)):
                    for j in range(len(result[0])):
                        finaltable[i][j+1]=result[i][j]
                optimalvalue="Does Not Exist"
                optimal_sol="Does Not Exist"
                my_dict={}
                my_dict['Initial Tableau']=initialtable
                my_dict['Final Tableau']=finaltable
                my_dict['Status of Solution']=statusofsol
                my_dict['Optimal Solution Vector']=optimal_sol
                my_dict['Optimal Value']=optimalvalue
                return my_dict
        elif cost!=0:
            statusofsol="Infeasible"
            finaltable=np.zeros((len(b),len(result[0])+1))
            for i in range(len(b)):
                finaltable[i][0]=b[i]
            for i in range(len(result)):
                for j in range(len(result[0])):
                    finaltable[i][j+1]=result[i][j]
            optimalvalue="Does Not Exist"
            optimal_sol="Does Not Exist"
            my_dict={}
            my_dict['Initial Tableau']=initialtable
            my_dict['Final Tableau']=finaltable
            my_dict['Status of Solution']=statusofsol
            my_dict['Optimal Solution Vector']=optimal_sol
            my_dict['Optimal Value']=optimalvalue
            return my_dict
        else:

            flag=False
            y=len(basicvarindex)
            for l in range(y):
                if basicvarindex[l]>n:
                    for j in range(n):
                        if (result[l][j]!=0):
                            flag=True
                            for u in range(len(result)):
                                if (u!=l):
                                    result[u,:]=result[u,:]-result[l,:]*result[u,j]/result[l,j]
                            basicvarindex[l]=j+1
                            reduced_cost=reduced_cost-result[l,:]*reduced_cost[j]/result[l,j]
                            result[l,:]=result[l,:]/result[l,j]
                            break
                    if flag==False:
                        continue
                else:
                    continue

            while(True):
                flag=True
                for l in range(len(basicvarindex)):
                    if basicvarindex[l]>n:
                        flag=False
                        result = np.delete(result, l, axis=0)
                        basicvarindex = np.delete(basicvarindex, l, axis=0)
                        b = np.delete(b, l, axis=0)
                        break
                if flag==True:
                    break
            reduced_cost_new=np.zeros(n)

            for i in range(0,n):
                reduced_cost_new[i]=reduced_cost[i]
                
            result = np.delete(result, np.s_[-(len(reduced_cost)-n):], axis=1)
            reduced_cost=reduced_cost_new

            x=np.zeros(n)
            for i in range(0,len(basicvarindex)):
                x[int(basicvarindex[i]-1)]=b[i]
            cost=-1* np.dot(x,modc)
            cb=np.zeros(len(basicvarindex))
            for i in range(len(basicvarindex)):
                cb[i]=modc[int(basicvarindex[i])-1]
            for i in range(0,n):
                if i+1 in basicvarindex:
                    reduced_cost[i]=0
                else:
                    reduced_cost[i]=modc[i]-np.dot(cb,result[:,i])

            numvar=len(result[0])
            numbasic=len(b)
            while True:
                result=np.round(result,6)
                reduced_cost=np.round(reduced_cost,6)
                cost=np.round(cost,6)
                b=np.round(b,6)
                if np.any(reduced_cost<0):
                    pivotcol = np.where(reduced_cost < 0)[0][0]
                    column_positive = np.any(result[:,pivotcol] > 0)
                    if column_positive:
                        mini=2**30-1
                        for j in range(0,len(b)):
                            if result[:,pivotcol][j]>0:
                                mini=min(mini,b[j]/result[:,pivotcol][j])
                        for j in range(0,len(b)):
                            if result[:,pivotcol][j]>0:
                                if b[j]/result[:,pivotcol][j]==mini:
                                    index1=j
                                    break
                        for i in range(0,numbasic):
                            if (i!=index1):
                                alpha=-1*(result[:,pivotcol][i]/result[:,pivotcol][index1])
                                result[i]=result[i]+(alpha)*result[index1]
                                b[i]=b[i]+alpha*b[index1]
                        alpha1=(-1*(reduced_cost[pivotcol])/result[:,pivotcol][index1])
                        reduced_cost=reduced_cost+alpha1*result[index1]
                        cost=cost+alpha1*b[index1]
                        b[index1]=b[index1]/result[:,pivotcol][index1]
                        result[index1]=result[index1]/result[:,pivotcol][index1]
                        basicvarindex[index1]=pivotcol+1
                    else:
                        statusofsol="Unbounded"
                        finaltable=np.zeros((len(b),len(result[0])+1))
                        for i in range(len(b)):
                            finaltable[i][0]=b[i]
                        for i in range(len(result)):
                            for j in range(len(result[0])):
                                finaltable[i][j+1]=result[i][j]
                        if (objective=="maximize"):
                            optimalvalue="inf"
                        else:
                            optimalvalue="-inf"
                        optimal_sol="Does Not Exist"
                        my_dict={}
                        my_dict['Initial Tableau']=initialtable
                        my_dict['Final Tableau']=finaltable
                        my_dict['Status of Solution']=statusofsol
                        my_dict['Optimal Solution Vector']=optimal_sol
                        my_dict['Optimal Value']=optimalvalue
                        return my_dict
                else:
                    statusofsol="Optimal"
                    finaltable=np.zeros((len(b),len(result[0])+1))
                    for i in range(len(b)):
                        finaltable[i][0]=b[i]
                    for i in range(len(result)):
                        for j in range(len(result[0])):
                            finaltable[i][j+1]=result[i][j]
                    if objective=='minimize':
                        optimalvalue=-1*cost
                    else:
                        optimalvalue=cost
                    final_bfs=np.zeros(len(A[0]))
                    for i in range(len(basicvarindex)):
                        if (basicvarindex[i]-1<len(A[0])):
                            final_bfs[int(basicvarindex[i]-1)]=b[i]
                    optimal_sol=final_bfs
                    my_dict={}
                    my_dict['Initial Tableau']=initialtable
                    my_dict['Final Tableau']=finaltable
                    my_dict['Status of Solution']=statusofsol
                    my_dict['Optimal Solution Vector']=optimal_sol
                    my_dict['Optimal Value']=optimalvalue
                    return my_dict
                    # break
print(simplex_algo("test.txt"))