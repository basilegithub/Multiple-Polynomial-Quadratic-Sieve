# This is the file containing the functions handling the relations list

from utils import invmod

def find_cycle(graph, init):
    if init[0] == 1: return DFS(graph, init)

    path_p1_to_1 = DFS(graph, [1, init[0]])

    if path_p1_to_1 is not None:
        path_p2_to_1 = DFS(graph, [1, init[1]])

        if path_p2_to_1 is not None:
            index1, index2 = len(path_p1_to_1)-1, len(path_p2_to_1)-1
            while path_p1_to_1[index1] == path_p2_to_1[index2]:
                index1 -= 1
                index2 -= 1
            return path_p2_to_1[:index2+2] + path_p1_to_1[:index1+1][::-1]

    return DFS(graph, init)

# Finds cycles, aka relations, in the graph of the large primes partial relations
def DFS(graph,init):
    path = [init[1]]
    stack = []
    neighbors = [init[1]]

    tmp = graph[init[1]]
    
    if init[0] in tmp:
        return path+[init[0]]
    
    neighbors += tmp
    
    stack.append(neighbors)
    
    # While there is a potential path that has not been explored
    while len(stack):
        while len(stack) and len(stack[-1]) == 1:
            path.pop()
            stack.pop()
        if not len(stack): break

        next_node = stack[-1][-1]
        stack[-1].pop()
        path.append(next_node)
        neighbors = [next_node]

        if init[0] in graph[next_node]:
            return path+[init[0]]

        parent = stack[-1][0]

        for i in range(len(graph[next_node])):
            if graph[next_node][i] != parent: neighbors.append(graph[next_node][i])
            else:
                neighbors += graph[next_node][i+1:]
                break

        stack.append(neighbors)
        
# Given a cycle, constructs the full relations
# The graph is always acyclic: when we find a new cycle, the partial relation we just found is not added
# Thus, one new partial relation can only create one cycle at most
def combine(relations,smooth,partial_relations,partial_smooth,path,value,n):
    smoo,relation = value, value*value-n
    for i in range(len(path)-1):
        firstp, secondp = min(path[i],path[i+1]), max(path[i],path[i+1])
        
        smoo = smoo*partial_smooth[firstp][secondp]%n
        relation *= partial_relations[firstp][secondp][0]
        
    for p in path:
        smoo = smoo*invmod(p,n)%n
        relation //= (p*p)
    smooth.append(smoo)
    relations.append(relation)
    
# Master function, keeps the partial_relations and possible_smooth lists sorted, find the matching large primes, create the full relations from large primes
def handle_possible_smooth(value, tmp_smooth, full_found, partial_found, relations, smooth_number, partial_relations, possible_smooth, graph, size_partials, connected_components, node_component, index_component, cycle_len, n):
    if tmp_smooth[0] == True:
        full_found += 1
        relations.append(value*value-n)
        smooth_number.append(value)

    elif tmp_smooth[0] == "large":
        if tmp_smooth[1] == tmp_smooth[2]:
            p = tmp_smooth[1]
            smoo = value*invmod(p,n)%n
            relation = (value*value-n)//(p*p)

            smooth_number.append(smoo)
            relations.append(relation)
            full_found += 1

        # if this is the first partial relation we see, there is none to combine it with
        elif not bool(partial_relations):
            pair = [value*value-n, tmp_smooth[1], tmp_smooth[2]]
            
            partial_relations[tmp_smooth[1]] = {}
            partial_relations[tmp_smooth[1]][tmp_smooth[2]] = pair
            
            possible_smooth[tmp_smooth[1]] = {}
            possible_smooth[tmp_smooth[1]][tmp_smooth[2]] = value
            
            graph[tmp_smooth[1]] = [tmp_smooth[2]]
            graph[tmp_smooth[2]] = [tmp_smooth[1]]

            connected_components[index_component] = set((tmp_smooth[1], tmp_smooth[2]))
            node_component[tmp_smooth[1]] = index_component
            node_component[tmp_smooth[2]] = index_component

            size_partials += 1
            index_component += 1
        else:
            pair = [value*value-n, tmp_smooth[1], tmp_smooth[2]]
            
            flag_small_prime = tmp_smooth[1] in graph

            if not flag_small_prime: # if we have never seen the small prime before
                graph[tmp_smooth[1]] = [tmp_smooth[2]]
                
                partial_relations[tmp_smooth[1]] = {}
                partial_relations[tmp_smooth[1]][tmp_smooth[2]] = pair

                possible_smooth[tmp_smooth[1]] = {}
                possible_smooth[tmp_smooth[1]][tmp_smooth[2]] = value
                
                flag_big_prime = tmp_smooth[2] in graph

                if flag_big_prime: # if we have seen the big prime before
                    graph[tmp_smooth[2]].append(tmp_smooth[1])

                    connected_components[node_component[tmp_smooth[2]]].add(tmp_smooth[1])
                    node_component[tmp_smooth[1]] = node_component[tmp_smooth[2]]

                else: # if we have not seen the big prime before
                    graph[tmp_smooth[2]] = [tmp_smooth[1]]

                    connected_components[index_component] = set((tmp_smooth[1], tmp_smooth[2]))
                    node_component[tmp_smooth[1]] = index_component
                    node_component[tmp_smooth[2]] = index_component

                    index_component += 1
               
            # If the smallest prime has already been seen
            else:
                flag_big_prime = tmp_smooth[2] in graph
                
                # If we have not seen the big prime before
                if not flag_big_prime:
                    graph[tmp_smooth[1]].append(tmp_smooth[2])
                    graph[tmp_smooth[2]] = [tmp_smooth[1]]

                    connected_components[node_component[tmp_smooth[1]]].add(tmp_smooth[2])
                    node_component[tmp_smooth[2]] = node_component[tmp_smooth[1]]
    
                    if tmp_smooth[1] not in partial_relations:
                        partial_relations[tmp_smooth[1]] = {}
                        possible_smooth[tmp_smooth[1]] = {}

                    partial_relations[tmp_smooth[1]][tmp_smooth[2]] = pair
                    possible_smooth[tmp_smooth[1]][tmp_smooth[2]] = value
                    size_partials += 1

                # If the largest prime has been seen, ie if both primes have already been seen
                else:
                    if node_component[tmp_smooth[1]] != node_component[tmp_smooth[2]]:
                        graph[tmp_smooth[1]].append(tmp_smooth[2])
                        graph[tmp_smooth[2]].append(tmp_smooth[1])

                        connected_components[node_component[tmp_smooth[1]]] = connected_components[node_component[tmp_smooth[1]]].union(connected_components[node_component[tmp_smooth[2]]])

                        for node in connected_components[node_component[tmp_smooth[2]]]:
                            node_component[node] = node_component[tmp_smooth[1]]

                        if tmp_smooth[1] not in partial_relations:
                            partial_relations[tmp_smooth[1]] = {}
                            possible_smooth[tmp_smooth[1]] = {}

                        partial_relations[tmp_smooth[1]][tmp_smooth[2]] = pair
                        possible_smooth[tmp_smooth[1]][tmp_smooth[2]] = value

                        size_partials += 1
                    else:
                        path_cycle = find_cycle(graph,[tmp_smooth[1], tmp_smooth[2]])
                        if len(path_cycle) < 11: cycle_len[len(path_cycle)-2] += 1
                        else: cycle_len[-1] += 1
                        combine(relations,smooth_number,partial_relations,possible_smooth,path_cycle,value,n)
                        partial_found += 1
                
    return relations, smooth_number, partial_relations, possible_smooth, full_found, partial_found, graph, size_partials, connected_components, node_component, index_component