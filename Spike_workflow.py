import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import Bio as Bio
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import sklearn
from sklearn.preprocessing import OneHotEncoder
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from plotly import subplots
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
from sklearn.manifold import MDS
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    #analysis_on_POI()
    graph_heatmap()
    #compare_to_wt()
    #make_hist()


def analysis_on_POI():
    matrix = {}
    count = 0
    for seq in SeqIO.parse('Whole_genome_with_B.1.1.7_patients_of_interest.fasta', format='fasta'):
        count += 1
        total_base_mismatches = []
        for seq_2 in SeqIO.parse('Whole_genome_with_B.1.1.7_patients_of_interest.fasta', format='fasta'):
            if seq_2.id != seq.id:
                print(seq.id, seq_2.id)
            count_of_mismatches = 0
            # spike is actually from 21563 to 25384 but 100bp upstream and downstream for errors
            for seq_bases, seq_2_bases in zip(seq.seq[21563:25384],
                                              seq_2.seq[21563:25384]):
                if seq_bases != seq_2_bases:
                    count_of_mismatches += 1
                else:
                    pass
            total_base_mismatches.append(count_of_mismatches)
        matrix[f'{seq.id}_{count}'] = total_base_mismatches
    df = pd.DataFrame.from_dict(matrix, orient='index')
    df.to_csv('Heatmap_data_for_Spike_.csv')


def graph_heatmap():
    df = pd.DataFrame(pd.read_csv('Heatmap_data_for_Spike_.csv'))
    labels = list(df.loc[:, 'ID'])
    df = df.drop(columns='ID')
    # h = sns.heatmap(df, linewidth=0.1, cmap=sns.color_palette('coolwarm', as_cmap=True),
    #               norm=norm, robust=True, fmt="f", xticklabels=labels, yticklabels=labels)

    h = sns.heatmap(df, vmax=200, vmin=0, center=100, linewidth=0.1,
                    yticklabels=labels, xticklabels=labels, cmap=sns.color_palette('rocket_r', as_cmap=True),
                    ax=sns.set(font_scale=0.35))
    plt.show()
    plt.savefig('Spike_DNA_heatmap.svg', format='svg')


def compare_to_wt():
    # todo getting many incomplete sequences can correct for it, but need to check in

    wt = 0
    for seq in SeqIO.parse('Whole_genome_with_B.1.1.7_patients_of_interest.fasta', format='fasta'):
        if seq.id == "B.1.1.7":
            wt = seq.seq
    wt = wt[21563:25384]  # todo ask nuno about this....
    matrix_of_bp = {}
    for wt_index, wt_element in enumerate(wt):
        count = 0
        for seq_ in SeqIO.parse('Whole_genome_with_B.1.1.7_patients_of_interest.fasta', format='fasta'):
            seq_trimmed = seq_.seq[21563:25384]
            seq_base = seq_trimmed[wt_index]
            if seq_base != wt_element:
                count += 1
            else:
                pass
        matrix_of_bp[wt_index] = count

    df = pd.DataFrame.from_dict(matrix_of_bp, orient='index')

    df.to_csv('Spike_mutation_compared_to_B.1.1.7.csv')


def make_hist():
    df_ = pd.read_csv('Spike_mutation_compared_to_B.1.1.7.csv')
    # drop last 500 positions due to maybe sequencing?????

    df_.columns = ['Base Position', 'Total # of Mutations compared to B.1.1.7']
    print(df_)
    hist_1 = px.histogram(df_, x='Base Position', y='Total # of Mutations compared to B.1.1.7',
                          labels={'x': 'Mutation', 'y': 'Instances of Mutation'})

    hist_1.write_html('Spike_mutations_to_B.1.1.7.html')

main()