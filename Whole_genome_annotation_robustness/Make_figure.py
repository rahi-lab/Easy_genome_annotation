import numpy as np
import matplotlib.pyplot as plt


font = {'fontname':'FreeSans'}

data = np.load('/mnt/DATA/Easy_genome_annotation/results_array.npy')  # shape (3, 100, 41)

mapped = data[0]  
low_score = data[1]
unmapped = data[2]


percent_indices = np.arange(mapped.shape[1]) / 2


fig, axes = plt.subplots(3, figsize=(10, 20), sharex=True)
fig.suptitle('Impact of Nucleotide Substitution on Read Mapping Performance', fontsize=32, **font)


plt.rc('xtick', labelsize=15) 
plt.rc('ytick', labelsize=17) 
titles = ['Mapped', 'Low MapQ Score', 'Not Found']
data_types = [mapped, low_score, unmapped]



for i, (ax, title, data_type) in enumerate(zip(axes, titles, data_types)):
    for iteration in range(data_type.shape[0]):
        ax.scatter(percent_indices, data_type[iteration], s = 10, alpha=0.03, c ='k')
        
    means = np.mean(data_type, axis=0) 
    standard_deviation_values  = np.std(data_type, axis=0, ddof=1)
    ax.errorbar(
        percent_indices, means, yerr=standard_deviation_values, fmt='o', capsize=4, markersize=1, elinewidth = 1, color='red', ecolor='red'         
    )
    ax.set_title(title,**font, fontsize=30)
    ax.set_ylabel('Read Count', **font, fontsize=22)
    


axes[-1].set_xlabel('Probability of Nucleotide Substitution (%)', **font, fontsize=22)  

plt.savefig("EGA_fig.svg", bbox_inches='tight')
plt.show()
