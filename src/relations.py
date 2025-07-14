# This is the file containing the functions handling the relations list

from utils import invmod

# Finds cycles, aka relations, in the graph of the large primes partial relations
def find_cycle(graph,init):
    path = [init[2]]
    next = []
    TMP = [init[2]]
    
    a,b = 0,len(graph)-1
    c = (a+b)>>1
    # Find where the largest prime is located in the graph
    while a <= b and graph[c][0] != init[2]:
        if graph[c][0] < init[2]: a = c+1
        else: b = c-1
        c = (a+b)>>1
    
    # Initialize the potential paths
    for i in range(1,len(graph[c])):
        if graph[c][i] == init[1]: return path+[init[1]]
        else: TMP.append(graph[c][i])
    
    next.append(TMP)
    
    # While there is a potential path that has not been explored
    while len(next):
        while len(next) and len(next[-1]) == 1:
            path.pop()
            next.pop()
        if not len(next): break
        goto = next[-1][-1]
        next[-1].pop()
        path.append(goto)
        TMP = [goto]
        if goto == 1: c = 0
        else:
            a,b = 0,len(graph)-1
            c = (a+b)>>1
            while a <= b and graph[c][0] != goto:
                if graph[c][0] < goto: a = c+1
                else: b = c-1
                c = (a+b)>>1
        for i in range(1,len(graph[c])):
            if graph[c][i] == init[1]: return path+[init[1]]
            elif graph[c][i] != next[-1][0]: TMP.append(graph[c][i])
        next.append(TMP)
        
# Given a cycle, constructs the full relations
# The graph is always acyclic: when we find a new cycle, the partial relation we just found is not added
# Thus, one new partial relation can only create one cycle at most
def combine(relations,smooth,partial_relations,partial_smooth,path,value,n):
    smoo,relation = value, value*value-n
    for i in range(len(path)-1):
        firstp,secondp = min(path[i],path[i+1]),max(path[i],path[i+1])
        if firstp == 1: a = 0
        else:
            a,b = 0,len(partial_relations)-1
            c = (a+b)>>1
            while a <= b:
                if partial_relations[c][1] < firstp: a = c+1
                else: b = c-1
                c = (a+b)>>1
        while a < len(partial_relations) and partial_relations[a][1] == firstp:
            if partial_relations[a][2] == secondp:
                smoo = smoo*partial_smooth[a]%n
                relation *= partial_relations[a][0]
                break
            else: a += 1
    for p in path:
        smoo = smoo*invmod(p,n)%n
        relation //= (p*p)
    smooth.append(smoo)
    relations.append(relation)
    
# Master function, keeps the partial_relations and possible_smooth lists sorted, find the matching large primes, create the full relations from large primes
def handle_possible_smooth(value, tmp_smooth, full_found, partial_found, relations, smooth_number, partial_relations, possible_smooth, graph, need_append, cycle_len, n):
    if tmp_smooth[0] == True:
        full_found += 1
        relations.append(value*value-n)
        smooth_number.append(value)
    elif tmp_smooth[0] == "large":
        #return relations, smooth_number, partial_relations, possible_smooth, need_append, full_found, partial_found
        if need_append:
            partial_relations.append([value*value-n,tmp_smooth[1], tmp_smooth[2]])
            possible_smooth.append(value)
            graph.append([tmp_smooth[1],tmp_smooth[2]])
            graph.append([tmp_smooth[2],tmp_smooth[1]])
            need_append = False
        else:
            A1,B1 = 0,len(graph)-1
            C1 = (A1+B1)>>1
            # Find in graph if smallest prime has already be seen
            while A1 <= B1 and graph[C1][0] != tmp_smooth[1]:
                if graph[C1][0] < tmp_smooth[1]: A1 = C1+1
                else: B1 = C1-1
                C1 = (A1+B1)>>1
            
            # If the smallest prime has not already been seen
            if C1 == -1  or graph[C1][0] != tmp_smooth[1]:
                graph.insert(A1,[tmp_smooth[1],tmp_smooth[2]]) # Insert smallest prime in graph
                # Find where in partial relations it is needed to insert the new partial relation
                A2,B2 = 0,len(partial_relations)-1
                C2 = (A2+B2)>>1
                while A2 <= B2:
                    if partial_relations[C2][1] < tmp_smooth[1]: A2 = C2+1
                    else: B2 = C2-1
                    C2 = (A2+B2)>>1
                partial_relations.insert(A2,[value*value-n,tmp_smooth[1],tmp_smooth[2]])
                possible_smooth.insert(A2,value)
                
                # Find if largest prime has already been seen
                A3,B3 = 0,len(graph)-1
                C3 = (A3+B3)>>1
                while A3 <= B3 and graph[C3][0] != tmp_smooth[2]:
                    if graph[C3][0] < tmp_smooth[2]: A3 = C3+1
                    else: B3 = C3-1
                    C3 = (A3+B3)>>1
                    
                # If it has not been seen, append it
                if C3 == -1 or graph[C3][0] != tmp_smooth[2]: graph.insert(A3,[tmp_smooth[2],tmp_smooth[1]])
                # Otherwise, add the smallest prime to the list of links of the largest prime
                else: graph[C3].append(tmp_smooth[1])
               
            # If the smallest prime has already been seen
            else:
                A2,B2 = 0,len(graph)-1
                C2 = (A2+B2)>>1
                # Find if the largest prime has already been seen
                while A2 <= B2 and graph[C2][0] != tmp_smooth[2]:
                    if graph[C2][0] < tmp_smooth[2]: A2 = C2+1
                    else: B2 = C2-1
                    C2 = (A2+B2)>>1
                
                # If it has not been seen
                if C2 == -1 or graph[C2][0] != tmp_smooth[2]:
                    graph[C1].append(tmp_smooth[2]) # Append the largest prime to the list of links of the smallest prime
                    graph.insert(A2,[tmp_smooth[2],tmp_smooth[1]]) # Insert the large prime in the graph
                    
                    # Find where to insert the partial relation
                    A3,B3 = 0,len(partial_relations)-1
                    C3 = (A3+B3)>>1
                    while A3 <= B3:
                        if partial_relations[C3][1] < tmp_smooth[1]: A3 = C3+1
                        else: B3 = C3-1
                        C3 = (A3+B3)>>1
                    partial_relations.insert(A3,[value*value-n,tmp_smooth[1],tmp_smooth[2]])
                    possible_smooth.insert(A3,value)

                # If the largest prime has been seen, ie if both primes have already been seen
                else:
                    path_cycle = find_cycle(graph,tmp_smooth)
                    if path_cycle == None: # If no cycle found
                        graph[C1].append(tmp_smooth[2])
                        graph[C2].append(tmp_smooth[1])
                        A3,B3 = 0,len(partial_relations)-1
                        C3 = (A3+B3)>>1
                        while A3 <= B3:
                            if partial_relations[C3][1] < tmp_smooth[1]: A3 = C3+1
                            else: B3 = C3-1
                            C3 = (A3+B3)>>1
                        partial_relations.insert(A3,[value*value-n,tmp_smooth[1],tmp_smooth[2]])
                        possible_smooth.insert(A3,value)
                    else:
                        if len(path_cycle) < 11: cycle_len[len(path_cycle)-2] += 1
                        else: cycle_len[-1] += 1
                        combine(relations,smooth_number,partial_relations,possible_smooth,path_cycle,value,n)
                        partial_found += 1
                
    return relations, smooth_number, partial_relations, possible_smooth, need_append, full_found, partial_found