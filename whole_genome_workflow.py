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
    # working_with_updated_alignemt()
    # analysis_on_POI()
    # graph_heatmap()
    #compare_to_wt()
    make_hist()

def alignment_name_change():
    alignment = AlignIO.read('whole_genome_with_B.1.1.7_alignment.aln', format='clustal')
    df = pd.DataFrame(pd.read_csv('Sequencing_metadata.csv'))
    list_of_files = []
    for file in alignment:  # iterate through all sequences in the alignment file
        file_id = file.id  # get the current ID
        if file_id == 'OW998408.1':
            file.id = 'B.1.1.7 Reference'
            list_of_files.append(file)
        else:
            print(file_id[-10:])
            room_and_patient = df.loc[df['NGI ID'] == file_id[-10:], ['Patient ID', 'Room ID', 'Date', 'Type']]
            patient = list(room_and_patient['Patient ID'])
            print(room_and_patient)
            room = list(room_and_patient['Room ID'])
            date_ = list(room_and_patient['Date'])
            type_ = list(room_and_patient['Type'])
            id_ = patient + room + date_ + type_
            print(id_)
            file.id = f'Patient_{id_[0]}_sample_{id_[3]}'
            print(file.id)
            list_of_files.append(file)
    SeqIO.write(list_of_files, 'Whole_genome_alignment_with_B.1.1.7.fasta', 'fasta')


def working_with_updated_alignemt():
    # generates new alignment file based on patients of interest
    list_of_pOI = []
    samples_ = ["Patient_1_", "Patient_2_", "Patient_3_",
                "Patient_4_", "Patient_6_", "Patient_7_", "Patient_9_",
                "Patient_12_", "Patient_13_", "Patient_14_", "Patient_15_", "Patient_16_",
                "Patient_18_", "Patient_20_", "Patient_27_", "Patient_34_", "B.1.1.7"]
    for seq in SeqIO.parse('Whole_genome_alignment_with_B.1.1.7.fasta', 'fasta'):
        for soi in samples_:
            print(seq.id)
            print(soi)
            if soi in str(seq.id):
                list_of_pOI.append(seq)
            else:
                pass
    SeqIO.write(list_of_pOI, 'Whole_genome_with_B.1.1.7_patients_of_interest.fasta', format='fasta')


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
            for seq_bases, seq_2_bases in zip(seq.seq[55:],
                                              seq_2.seq[55:]):  # todo check with Nuno if primers were trimmed
                if seq_bases != seq_2_bases:
                    count_of_mismatches += 1
                else:
                    pass
            total_base_mismatches.append(count_of_mismatches)
        matrix[f'{seq.id}_{count}'] = total_base_mismatches
    df = pd.DataFrame.from_dict(matrix, orient='index')
    df.to_csv('Heatmap_data_for_whole_genome_.csv')


def compare_to_wt():
    # todo getting many incomplete sequences can correct for it, but need to check in
    wt = 0
    for seq in SeqIO.parse('Whole_genome_with_B.1.1.7_patients_of_interest.fasta', format='fasta'):
        if seq.id == "B.1.1.7":
            wt = seq.seq
    wt = wt[55:]  # todo ask nuno about this....
    matrix_of_bp = {}
    for wt_index, wt_element in enumerate(wt):
        count = 0
        for seq_ in SeqIO.parse('Whole_genome_with_B.1.1.7_patients_of_interest.fasta', format='fasta'):
            seq_trimmed = seq_.seq[55:]
            seq_base = seq_trimmed[wt_index]
            if seq_base != wt_element:
                count += 1
            else:
                pass
        matrix_of_bp[wt_index] = count

    df = pd.DataFrame.from_dict(matrix_of_bp, orient='index')
    df_1 = df.iloc[:15000, :]
    df_2 = df.iloc[15001:, :]

    df_1.to_csv('Whole_genome_mutation_0.1_compared_to_B.1.1.7.csv')
    df_2.to_csv('Whole_genome_mutation_0.2_compared_to_B.1.1.7.csv')


def make_hist():
    df_1 = pd.read_csv('Whole_genome_mutation_0.1_compared_to_B.1.1.7.csv')
    df_2 = pd.read_csv('Whole_genome_mutation_0.2_compared_to_B.1.1.7.csv')
    list_of_df = [df_1, df_2]
    df_ = pd.concat(list_of_df)
    # drop last 500 positions due to maybe sequencing?????

    df_.columns = ['Base Position', 'Total # of Mutations compared to B.1.1.7']
    df_ = df_.iloc[:-55, :]
    print(df_)
    hist_1 = px.histogram(df_, x='Base Position', y='Total # of Mutations compared to B.1.1.7',
                          labels={'x': 'Mutation', 'y': 'Instances of Mutation'})

    hist_1.write_html('Whole_genome_mutations_to_B.1.1.7.html')



def graph_heatmap():
    df = pd.DataFrame(pd.read_csv('Heatmap_data_for_whole_genome_.csv'))
    labels = list(df.loc[:, 'ID'])
    df = df.drop(columns='ID')
    # h = sns.heatmap(df, linewidth=0.1, cmap=sns.color_palette('coolwarm', as_cmap=True),
    #               norm=norm, robust=True, fmt="f", xticklabels=labels, yticklabels=labels)

    h = sns.heatmap(df, vmax=400, vmin=0, center=200, linewidth=0.1,
                    yticklabels=labels, xticklabels=labels, cmap=sns.color_palette('rocket_r', as_cmap=True),
                    ax=sns.set(font_scale=0.35))
    plt.show()
    plt.savefig('Whole_genome_DNA_heatmap.svg', format='svg')


main()
