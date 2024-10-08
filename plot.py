import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



ens = 50    #number of ensemble
N = 16
T_c = np.sqrt(N-1)

data_agr_I = pd.read_csv(f"Coevolutionary/Coev_16_agreement_1.csv")
data_agr_II = pd.read_csv(f"Coevolutionary/Coev_{N}_agreement_-1.csv")

data_agr_I["scale_temp"] = (data_agr_I["Temperature"] - T_c) / T_c
data_agr_II["scale_temp"] = (data_agr_II["Temperature"] - T_c) / T_c

T = data_agr_I["scale_temp"]
E = data_agr_I["Energy"]
Q_agr_I = data_agr_I["node link corr"]
Q_agr_II = data_agr_II["node link corr"]

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


plt.savefig(f'Coevolutionary/Coev_{N}_ensemble_{ens}_cpp.png', dpi=300)
plt.close()