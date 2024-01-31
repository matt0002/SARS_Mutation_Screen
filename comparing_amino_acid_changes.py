import pandas as pd
import plotly.express as px


def format():
    df_in = pd.DataFrame(pd.read_csv('Whole_genome_to_amino_acid.csv'))

    col_of_interest = 'numMutsXonly_list'
    col_of_interest_1 = 'UniqueMutList'
    col_of_interest_2 = 'ExistingMutList'

    dict_updated = {}
    samples_of_interest = ["Patient_1", "Patient_2", "Patient_3",
                           "Patient_4", "Patient_6", "Patient_7", "Patient_9",
                           "Patient_12", "Patient_13", "Patient_14", "Pattient_15", "Patient_16",
                           "Patient_18", "Patient_20", "Patient_27", "Patient_34", "B.1.1.7"]

    all_patients = list(df_in.loc[:, 'Patient ID'])

    for samples in all_patients:
        for soi in samples_of_interest:
            if soi in samples:
                room_ = list(df_in.loc[df_in['Patient ID'] == samples, 'Room ID'])
                date_ = list(df_in.loc[df_in['Patient ID'] == samples, 'Date'])
                muts = list(df_in.loc[df_in['Patient ID'] == samples, col_of_interest])
                muts_1 = list(df_in.loc[df_in['Patient ID'] == samples, col_of_interest_1])
                muts_2 = list(df_in.loc[df_in['Patient ID'] == samples, col_of_interest_2])

                all_muts = muts + muts_1 + muts_2

                name_ = f'{samples}_room_{room_[0]}_{date_[0]}'
                dict_updated[name_] = all_muts

    print(dict_updated)

    dict_unique_index = {}
    for key_ in dict_updated.keys():
        try:
            vals = dict_updated[key_]
            for val in vals:
                new_val = list(val.split(','))
                dict_unique_index[key_] = new_val
        except AttributeError:
            pass

    df_out = pd.DataFrame.from_dict(dict_unique_index, orient='index')
    df_out = df_out.transpose()

    df_out.to_csv('Amino_acid_whole_genome_unique_mutations.csv', index=False)



def plotting_site_occurrences():
    df = pd.DataFrame(pd.read_csv('Amino_acid_whole_genome_unique_mutations.csv'))
    # idea to break up mutation scheme by types of mutations and where they are

    # first get all sites

    all_mutation_sites = []
    for col in df.columns:
        try:
            all_muts = list(df.loc[:, col])
            for mut in all_muts:
                if type(mut) == float:
                    pass
                else:
                    if '(' in mut:
                        mut = mut[1:]
                    elif ')' in mut:
                        mut = mut[:-1]
                new_split = mut.split('_')
                if new_split[0] in all_mutation_sites:
                    pass
                else:
                    all_mutation_sites.append(new_split[0])
        except AttributeError:
            pass

    dict_of_mut_occurances = {}

    for unique_mut in all_mutation_sites:
        try:
            count = 0
            for col in df.columns:
                sites_ = list(df.iloc[:, col])
                for mut_by_sample in sites_:
                    if type(mut_by_sample) == float:
                        pass
                    else:
                        if unique_mut in mut_by_sample:
                            count += 1
            dict_of_mut_occurances[unique_mut] = count
        except AttributeError:
            pass

    print(dict_of_mut_occurances)
    hist = px.histogram(x=dict_of_mut_occurances.keys(), y=dict_of_mut_occurances.values(), )
    hist.write_html('hist_of_whole_genome_mutation_occurrences.html')


def mutations_in_site():
    # idea is to break up mutations by their protein
    df = pd.DataFrame(pd.read_csv('Amino_acid_whole_genome_unique_mutations.csv'))
    all_mutation_sites = []
    for col in df.columns:
        try:
            all_muts = list(df.loc[:, col])
            for mut in all_muts:
                if type(mut) == float:
                    pass
                else:
                    if '(' in mut:
                        mut = mut[1:]
                    elif ')' in mut:
                        mut = mut[:-1]
                new_split = mut.split('_')
                if new_split[0] in all_mutation_sites:
                    pass
                else:
                    all_mutation_sites.append(new_split[0])
        except AttributeError:
            pass

    dict_of_site_muts = {}

    for unique_mutations in all_mutation_sites:
        list_of_muts = []
        for col in df.columns:
            sites_ = list(df.loc[:, col])
            for mutation in sites_:
                mutation = str(mutation)
                print(mutation)
                if unique_mutations in mutation:
                    split_ = mutation.split('_')
                    list_of_muts.append(split_[1])

        dict_of_site_muts[unique_mutations] = list_of_muts

    for protein in dict_of_site_muts.keys():
        mutations_in_protein = dict_of_site_muts[protein]
        dict_counts = {}
        for individual_mutation in mutations_in_protein:
            count = 0
            for compare_mutations in mutations_in_protein:
                if individual_mutation == compare_mutations:
                    count += 1
            dict_counts[individual_mutation] = count
        print(dict_counts)
        if len(dict_counts.keys()) == 0:
            pass
        else:
            hist = px.histogram(x=dict_counts.keys(), y=dict_counts.values(), labels={'x': 'Mutations',
                                                                                      'y': 'Occurances seen across Patients'})
            hist.write_html(f'Mutation_occurances_in_{protein}.html')

format()
mutations_in_site()
