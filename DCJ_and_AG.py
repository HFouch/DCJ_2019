import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import copy

genomeA = [[ '.',1,4,-7,2,-6,-5,3,'.']]
genomeB = [['.', 1,2,3,4,5,6,7,'.']]

def add_anchors(genome):
    #NEED DIFFERENT ANCHORING STRATEGY IF IT IS TO WORK FOR MULTICHROMOSOMAL GENOMES.
    #Ask Riaan for input -he mentioned something in a lab meeting

    genes = genome[0][1:-1]
    anchored_genome = ['.']
    anchored_genome.append(0)
    for gene in genes:
        anchored_genome.append(gene)
    anchored_genome.append(len(genes)+1)
    anchored_genome.append('.')
    return anchored_genome

def create_gene_extremities(anchored_genome):
    gene_extremities = []
    for marker in anchored_genome:
        #if marker == '.'or marker == '^':
         #   gene_extremities.append(marker)
        if marker == '.':
            gene_extremities.append('-')


        elif int(marker) >= 0:
            marker_str = str(marker)
            gene_extremities.append(marker_str+'t')
            gene_extremities.append(marker_str+'h')
        else:
            marker_str = str(abs(marker))
            gene_extremities.append(marker_str + 'h')
            gene_extremities.append(marker_str + 't')

    return gene_extremities


def create_adjacencies(gene_extremities):
    adjacencies = []
    counter = 0
    while counter < len(gene_extremities):
        adjacencies.append((gene_extremities[counter], gene_extremities[counter + 1]))
        counter += 2
        '''

        if gene_extremities[counter] == gene_extremities[0]:
            adjacencies.append((gene_extremities[counter]))
            counter+=1
        elif gene_extremities[counter] == gene_extremities[-1]:
            adjacencies.append((gene_extremities[counter]))
            counter +=1

        else:
            adjacencies.append((gene_extremities[counter], gene_extremities[counter+1]))
            counter+=2

        '''
    return adjacencies

anchored_genomeA = (add_anchors(genomeA))
gene_extremitiesA = create_gene_extremities(anchored_genomeA)
adjacenciesA = create_adjacencies(gene_extremitiesA)
print('adj a: ', adjacenciesA)

anchored_genomeB = (add_anchors(genomeB))
gene_extremitiesB = create_gene_extremities(anchored_genomeB)
adjacenciesB = create_adjacencies(gene_extremitiesB)
print()
print('adj b: ', adjacenciesB)

def create_adjacency_table(adjacencies):
    adj_table_first = []
    adj_table_second =[]
    adj_table =[]
    i=0
    for adj in adjacencies:

        adj_table_first.append(adj[0])
        adj_table_second.append(adj[1])

        '''

        if len(adj) ==1:

            adj_table_first.append(adj)
            adj_table_second.append('-')
        else:

            adj_table_first.append(adj[0])
            adj_table_second.append(adj[1])
            
        '''

    adj_table.append(adj_table_first)
    adj_table.append(adj_table_second)
    return adj_table

adjacency_tableA = create_adjacency_table(adjacenciesA)
print(adjacency_tableA)
adjacency_tableB = create_adjacency_table(adjacenciesB)
print(adjacency_tableB)

#This is correct
print('Adjacency Table A')
print(adjacency_tableA[0])
print(adjacency_tableA[1])

def sort_one_level_down(adjacenciesB, adjacenciesA, adjacency_tableA):
    level_operations = []
    level_adj_tables = []
    adjacency_tableA_copy = copy.deepcopy(adjacency_tableA)

    for adj in adjacenciesB:
        print('###')
        print('adj: ', adj)
        intermediate_table = []
        print(intermediate_table)
        print(adjacency_tableA_copy[0])
        print(adjacency_tableA_copy[1])
        print()

        intermediate_table = copy.deepcopy(adjacency_tableA_copy)
        print(intermediate_table[0])
        print(intermediate_table[1])
        print()
        #If it is not a telomeric adjacency do:
        if adj[0] != '-' and adj[1] != '-':

            p = adj[0]
            q = adj[1]

            if p in intermediate_table[0]:
                p_index_in_A = intermediate_table[0].index(p)
            elif p in intermediate_table[1]:
                p_index_in_A = intermediate_table[1].index(p)
            u = ((intermediate_table[0][p_index_in_A], intermediate_table[1][p_index_in_A]))

            if q in intermediate_table[0]:
                q_index_in_A = intermediate_table[0].index(q)
            elif q in intermediate_table[1]:
                q_index_in_A = intermediate_table[1].index(q)
            v = ((intermediate_table[0][q_index_in_A], intermediate_table[1][q_index_in_A]))


            #NTS
            #still need to account for the extremities in the adjacency being switched around in intermediates

            #if the adjacency in B does not exist in A do:
            if u != v:

                if u[0] == p:
                    intermediate_table[0][p_index_in_A] = u[0]
                    intermediate_table[0][q_index_in_A] = u[1]

                elif u[1] == p:
                    intermediate_table[0][p_index_in_A] = u[1]
                    intermediate_table[0][q_index_in_A] = u[0]

                if v[0] == q:
                    intermediate_table[1][p_index_in_A] = v[0]
                    intermediate_table[1][q_index_in_A] = v[1]
                elif v[1] == q:
                    intermediate_table[1][p_index_in_A] = v[1]
                    intermediate_table[1][q_index_in_A] = v[0]

                operation = ((u, v), ((intermediate_table[0][p_index_in_A], intermediate_table[1][p_index_in_A]),
                                      (intermediate_table[0][q_index_in_A], intermediate_table[1][q_index_in_A])))


                print('operation: ',operation)
                print(intermediate_table[0])
                print(intermediate_table[1])
                print()

                level_adj_tables.append(intermediate_table)
                level_operations.append(operation)
                intermediate_table = []

        #If adj is a telomeric adjacency do:
        else:
            if adj[0] != '-':
                p = adj[0]
            else:
                p = adj[1]


            if p in adjacency_tableA[0]:
                p_index_in_A = adjacency_tableA[0].index(p)
            elif p in adjacency_tableA[1]:
                p_index_in_A = adjacency_tableA[1].index(p)
            u = ((adjacency_tableA[0][p_index_in_A], adjacency_tableA[1][p_index_in_A]))

            if u[0] == p:
                if u[1] != '-':
                    adjacency_tableA[0][p_index_in_A] = u[0]
                    adjacency_tableA[1][p_index_in_A] = '-'
                    adjacency_tableA[0].append('-')
                    adjacency_tableA[1].append(u[1])

                    operation = ((u), ((intermediate_table[0][p_index_in_A], intermediate_table[1][p_index_in_A]),
                                          (intermediate_table[0][q_index_in_A], intermediate_table[1][q_index_in_A])))

                    level_adj_tables.append(intermediate_table)
                    level_operations.append(operation)

            elif u[1] == p:
                if u[0] != '-':
                    adjacency_tableA[0][p_index_in_A] = '-'
                    adjacency_tableA[1][p_index_in_A] = u[1]
                    adjacency_tableA[0].append(u[0])
                    adjacency_tableA[1].append('-')

                    operation = ((u), ((intermediate_table[0][p_index_in_A], intermediate_table[1][p_index_in_A]),
                                       (intermediate_table[0][q_index_in_A], intermediate_table[1][q_index_in_A])))
                    level_adj_tables.append(intermediate_table)
                    level_operations.append(operation)




    return level_operations, level_adj_tables

operations_per_level = sort_one_level_down(adjacenciesB,adjacenciesA, adjacency_tableA)
print()
print('level operations')
print(operations_per_level[0])
print()
print('adj tables')
print(operations_per_level[1])
print()

print(adjacency_tableA[0])
print(adjacency_tableA[1])
print()
for i in range(0, len(operations_per_level[0])):
    print(operations_per_level[0][i])
    print(operations_per_level[1][i][0])
    print(operations_per_level[1][i][1])
    print()


def greedy_DCJ_sorting(adjacenciesB, adjacenciesA, adjacency_tableA):


    all_intermediates = []
    all_intermediates.append(adjacenciesA)



    for adj in adjacenciesB:



        if adj[0] != '-' and adj[1] != '-':

            p = adj[0]
            q = adj[1]

            if p in adjacency_tableA[0]:
                p_index_in_A = adjacency_tableA[0].index(p)
            elif p in adjacency_tableA[1]:
                p_index_in_A = adjacency_tableA[1].index(p)
            u = ((adjacency_tableA[0][p_index_in_A], adjacency_tableA[1][p_index_in_A]))


            if q in adjacency_tableA[0]:
                q_index_in_A = adjacency_tableA[0].index(q)
            elif q in adjacency_tableA[1]:
                q_index_in_A = adjacency_tableA[1].index(q)
            v = ((adjacency_tableA[0][q_index_in_A], adjacency_tableA[1][q_index_in_A]))


            if u != v:


                if u[0] == p:
                    adjacency_tableA[0][p_index_in_A] = u[0]
                    adjacency_tableA[0][q_index_in_A] = u[1]

                elif u[1] == p:
                    adjacency_tableA[0][p_index_in_A] = u[1]
                    adjacency_tableA[0][q_index_in_A] = u[0]

                if v[0] == q:
                    adjacency_tableA[1][p_index_in_A] = v[0]
                    adjacency_tableA[1][q_index_in_A] = v[1]
                elif v[1] == q:
                    adjacency_tableA[1][p_index_in_A] = v[1]
                    adjacency_tableA[1][q_index_in_A] = v[0]



                operation = ((u, v), ((adjacency_tableA[0][p_index_in_A], adjacency_tableA[1][p_index_in_A]),(adjacency_tableA[0][q_index_in_A], adjacency_tableA[1][q_index_in_A])))
                print()
                print('OPERATION: ', operation)
                print()

                intermediate_adjacenciesA = []
                for i in range(0, len(adjacency_tableA[0])):
                    intermediate_adjacenciesA.append((adjacency_tableA[0][i], adjacency_tableA[1][i]))

                all_intermediates.append(intermediate_adjacenciesA)
                print('inter adjs: ', intermediate_adjacenciesA)

        #else we are dealing with a telomere
        else:
            if adj[0] != '-':
                p = adj[0]
            else:
                p = adj[1]

            #need an if that ensures u is an adjacency

            if p in adjacency_tableA[0]:
                p_index_in_A = adjacency_tableA[0].index(p)
            elif p in adjacency_tableA[1]:
                p_index_in_A = adjacency_tableA[1].index(p)
            u = ((adjacency_tableA[0][p_index_in_A], adjacency_tableA[1][p_index_in_A]))

            if u[0] ==p:
                if u[1] !='-':
                    adjacency_tableA[0][p_index_in_A] = u[0]
                    adjacency_tableA[1][p_index_in_A] = '-'
                    adjacency_tableA[0].append('-')
                    adjacency_tableA[1].append(u[1])
            elif u[1] == p:
                if u[0] != '-':
                    adjacency_tableA[0][p_index_in_A]='-'
                    adjacency_tableA[1][p_index_in_A]=u[1]
                    adjacency_tableA[0].append(u[0])
                    adjacency_tableA[1].append('-')








    if len(all_intermediates) > 1:
        if adj == adjacenciesB[-1]:
            all_intermediates.pop()
            all_intermediates.append(adjacenciesB)
    #if len(all_intermediates) == 0:
     #   print('LEN == 0!')
      #  all_intermediates.append(adjacenciesB)
    return adjacency_tableA, all_intermediates


DCJ_sort = greedy_DCJ_sorting(adjacenciesB, adjacenciesA, adjacency_tableA)
sorted_adjacency_tableA = DCJ_sort[0]
print()
print('sorted table: ', sorted_adjacency_tableA)
print()
intermediary_adj_lists = DCJ_sort[1]
print(intermediary_adj_lists)
print()
print()
print('AdjacenciesB then AdjacenciesA final')
print(adjacenciesB)
print(intermediary_adj_lists[-1])
print()
print()



def add_genome_id_to_adjacencies(adjacenciesA, adjacenciesB):
    adjacenciesA_id = []
    adjacenciesB_id = []
    for (x, y) in adjacenciesA:
        new_adj = ('A'+x, 'A'+y)
        adjacenciesA_id.append(new_adj)

    for (x, y) in adjacenciesB:
        new_adj = ('B'+x, 'B'+y)
        adjacenciesB_id.append(new_adj)

    return adjacenciesA_id, adjacenciesB_id

id_adjacencies = add_genome_id_to_adjacencies(adjacenciesA,adjacenciesB)
adjacenciesA_id = id_adjacencies[0]
adjacenciesB_id = id_adjacencies[1]


def build_adj_graph(adjacenciesA, adjacenciesB, fig_name='adj_graph.png'):
    id_adjacencies = add_genome_id_to_adjacencies(adjacenciesA, adjacenciesB)
    adjacenciesA_id = id_adjacencies[0]
    adjacenciesB_id = id_adjacencies[1]
    print()
    print('building graph...')
    print('adjB: ', adjacenciesB_id)
    print('adjA: ', adjacenciesA_id)
    print()
    adj_graph = nx.Graph()
    for adj in adjacenciesA_id:
        adj_graph.add_node(adj, bipartite=0)
    for adj in adjacenciesB_id:
        adj_graph.add_node(adj, bipartite=1)

    edge_list = []


    for (x, y) in adjacenciesA:
        print(('x,y: ', (x, y)))
        for (a, b) in adjacenciesB:
            if x == a or x == b:
                b1 = (a, b)
                #edge_list.append((('A'+x, 'A'+y ),('B'+a, 'B'+b)))
                edge_list.append((adjacenciesA_id[adjacenciesA.index((x, y))], adjacenciesB_id[adjacenciesB.index((a,b))]))
                break
        for (c, d) in adjacenciesB:
            if y == c or y==d:
                b2 = (c, d)
                edge_list.append(
                    (adjacenciesA_id[adjacenciesA.index((x, y))], adjacenciesB_id[adjacenciesB.index((c,d))]))
                #edge_list.append((('A' + x, 'A' + y), ('B' + c, 'B' + d)))
                break



    connecting_edges = []

    for adj in adjacenciesB_id:
        connecting_edges.append((adjacenciesA_id[0], adj))
    for adj in adjacenciesA_id:
        connecting_edges.append((adj, adjacenciesB_id[0]))


    invisible_edges = []
    for edge in connecting_edges:
        if edge not in edge_list:
            invisible_edges.append(edge)

    '''
    connecting_edges = []
    connecting_edges.append((adjacenciesB_id[0], adjacenciesA_id[0]))
    connecting_edges.append((adjacenciesA_id[len(adjacenciesA_id)-1], adjacenciesB_id[len(adjacenciesB)-1]))
    for i in range(0, len(adjacenciesA_id)):
        if i == 0:
            connecting_edges.append((adjacenciesA_id[i], adjacenciesB_id[i+1]))
        elif i > len(adjacenciesB_id)-2:
            connecting_edges.append((adjacenciesA_id[i], adjacenciesB_id[i-1]))
        else:
            connecting_edges.append((adjacenciesA_id[i], adjacenciesB_id[i + 1]))
            connecting_edges.append((adjacenciesA_id[i], adjacenciesB_id[i - 1]))
    print('connecting edges: ')
    print(connecting_edges)
    invisible_edges = []
    for edge in connecting_edges:
        if edge not in edge_list:
            invisible_edges.append(edge)
    '''
    print()
    print('Invisible edges')
    print(invisible_edges)
    print()
    print('Edge list')
    print(edge_list)


    adj_graph.add_edges_from(edge_list)
    adj_graph.add_edges_from(invisible_edges)

    l, r = nx.bipartite.sets(adj_graph)
    pos = {}
    pos.update((node, (index, 1)) for index, node in enumerate(l))
    pos.update((node, (index, 2)) for index, node in enumerate(r))

    edge_colors = ['white' if edge in invisible_edges else 'black' for edge in adj_graph.edges()]
    nx.draw_networkx(adj_graph, pos=pos, with_labels=True, node_color='black', node_size=5, edge_color=edge_colors)
    #plt.show()
    #plt.savefig('save.png')
    #figure = plt.gcf()  # get current figure
    #figure.set_size_inches(8, 6)
    # when saving, specify the DPI
    plt.show()
    #nx.draw_networkx(adj_graph, pos=pos, with_labels=True, node_color='black', node_size=5, edge_color=edge_colors)
    #plt.savefig(fig_name, dpi=100,pad_inches=5)




def create_adj_graphs_for_intermediates(adjacenciesB, intermediary_adjacenciesA):

    for intermediate in intermediary_adjacenciesA:
        '''
        add_adjacencies_id = add_genome_id_to_adjacencies(intermediate, adjacenciesB)
        adjacenciesA_id = add_adjacencies_id[0]
        adjacenciesB_id = add_adjacencies_id[1]
        print()
        print('first b then a')
        print(adjacenciesB_id)
        print(adjacenciesA_id)
        print()
        '''
        if intermediate == intermediary_adjacenciesA[-1]:
            print('last intermediate')
            print(intermediate)
            print()
            fig_name = 'adjacency_graph_sorted_genome.png'
            build_adj_graph(intermediate, adjacenciesB, fig_name=fig_name)
        else:

            fig_name = 'adjacency_graph_intermediate' + str(intermediary_adjacenciesA.index(intermediate)+1)+'.png'
            build_adj_graph(intermediate, adjacenciesB, fig_name=fig_name)

#create_adj_graphs = create_adj_graphs_for_intermediates(adjacenciesB, intermediary_adj_lists)
'''
basically what you want to do next is make id'ed versions for b and each a and then make and save the adj graph for each of them 
def create_adj_graphs_for_intermediates(adjacenciesB_id, intermediary_adjacenciesA):
    for element in intermediary_adj_lists:
        adjacenciesA_id =
 '''

#print()
#for intermediate in intermediary_adj_lists:
#    print(intermediate)
#    print()
