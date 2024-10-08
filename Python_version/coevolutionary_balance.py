import numpy as np
import matplotlib.pyplot as plt
import random
import math
from numba import njit


@njit
def initial_random_matrix(N):
    """initial random adjacency matrix btx -1 and 1

    Args:
        N (int): matrix dimention (number of nodes)

    Returns:
        2d_array: adjacency matrix
    """
    A = np.zeros((N,N))
    for i in range(N):
        for j in range(i):
            A[i,j] = random.choice([-1,1])
            A[j,i] = A[i,j]
    return A

@njit
def initial_random_nodes(N):
    """initial random nodes attributes btw -1 and 1

    Args:
        N (int): number of nodes

    Returns:
        1d_array: nodes array
    """
    return random.choice([-1, 1], N)


def initial_agreement(N, agr=1):
    """inital adjacency matrix and nodes array from agreement type (default agreement is 1) 

    Args:
        N (int): number of nodes
        agr (int, optional): agreement. Defaults to 1.

    Returns:
        tuple: (adjacency matrix(N*N) , nodes array(N))
    """
    nodes = np.ones(N)
    A = np.ones((N,N)) - np.eye(N)
    if agr == -1:
        nodes *= -1
    return A , nodes



@njit
def check_flip(A, nodes, T):
    """check flip the node or link with uniform probability based on metropolis algorithm

    Args:
        A (2d_array): adajency matrix
        nodes (1d_array): nodes array
        T (float): temperature

    Returns:
        tuple: (adjacency matrix(N*N) , nodes array(N))
    """
    N = len(nodes)
    P = N/(N + N*(N - 1)/2)

    #link choose
    if random.random() > P:
        row , col = random.randint(0,N-1), random.randint(0,N-1)
        A_new = np.copy(A)
        A_new[row,col] *= -1
        A_new[col,row] *= -1
        dE = nodes[row]*nodes[col]*(- A_new[row,col] + A[row,col])
        
        if  dE <= 0:
            return A_new, nodes

        elif random.random() < math.exp(-dE/T):
            return A_new, nodes
        
        else:
            return A, nodes
    
    #node choose
    else:
        index = random.randint(0,N-1)
        nodes_new = np.copy(nodes)
        nodes_new[index] *= -1
        dE = (- nodes_new[index] + nodes[index]) * (np.sum(A[index,:]*nodes))

        if  dE <= 0:
            return A, nodes_new

        elif random.random() < math.exp(-dE/T):
            return A, nodes_new
            
        else:
            return A, nodes
            
@njit
def node_link_corr(A, nodes):
    """calculate node link corrolation

    Args:
        A (2d_array): adajency matrix
        nodes (1d_array): nodes array

    Returns:
        float: node link corr
    """
    return (nodes @ A)/len(nodes)

@njit
def node_link_node_corr(A, nodes):
    """calculate  the energy of network per number of triplets

    Args:
        A (2d_array): adajency matrix
        nodes (1d_array): nodes array

    Returns:
        float: meean (node link nodes) 
    """
    return -(nodes @ A @ nodes.T)/(len(nodes)**2)

def main(N , ens, T):

    Q_agr_I = []      #list of node link corr for initial agreement = 1
    Q_agr_II = []      #list of node link corr
    E = []      #list of node link node corr


    for t in T:
        Q_ens_I, Q_ens_II, E_ens = [], [], []
        for e in range(ens):
            A_I , nodes_I = initial_agreement(N, agr=1)
            A_II , nodes_II = initial_agreement(N, agr=-1)
            for _ in range(N**3):
                A_I , nodes_I = check_flip(A_I , nodes_I, t)
                A_II , nodes_II = check_flip(A_II , nodes_II, t)
            Q_ens_I.append(node_link_corr(A_I, nodes_I))
            Q_ens_II.append(node_link_corr(A_II, nodes_II))
            E_ens.append(node_link_node_corr(A_I, nodes_I))
        Q_agr_I.append(np.mean(Q_ens_I))
        Q_agr_II.append(np.mean(Q_ens_II))
        E.append(np.mean(E_ens))


    T_c = np.sqrt(N-1)
    T = (T - T_c) / T_c

    return Q_agr_I, Q_agr_II , E , T



#main program
if __name__ == "__main__":

    ens = 50    #number of ensemble
    N = 16
    T = np.linspace(0.01, 10 , 40)     #array of temperature
    Q_agr_I, Q_agr_II , E, T = main(N, ens, T)

    fig, axs = plt.subplots(1, 2, figsize= (16,8))
    fig.tight_layout(pad=4.0)

    for ax in axs.flatten():
        ax.tick_params(color='#708090', labelcolor='#708090')
        for spine in ax.spines.values():
            spine.set_edgecolor('#708090')

    fig.suptitle(f'Number of nodes is {N} and number of ensemble {ens}', fontsize=16)

    axs[0].plot(T, Q_agr_I, marker='o', color='#661D98',ls='None')
    axs[0].plot(T, Q_agr_II, marker='o', color='#661D98',ls='None')
    axs[0].set(xlabel= 'Temperature', ylabel= r'$<s_{i}\sigma_{ij}>$')
    axs[0].set_title('The node-link correlation')


    axs[1].plot(T, E, marker='o', color='#2C6FFE',ls='None')
    axs[1].set(xlabel= 'Temperature ', ylabel= r'$-<s_{i}\sigma_{ij}s_j>$')
    axs[1].set_title('The energy of network per number of triplets')


    plt.savefig(f'Coevolutionary/Coev_{N}_ensemble_{ens}.png', dpi=300)
    plt.close()
