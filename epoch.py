"""
Extracting Drugs from Clinical Trials
Author: Prabu Siva

Each task has it's own function:
Task 1: match_drug_names() returns ./output/matched_drug_names.json
Task 2: match_usan_codes() returns ./output/drugs_usan.json
Task 3: counts_of_trials_by_usan_class() returns ./output/trials_by_usan.json
Task 4: agg_counts_of_usan_pairs() return ./output/counts_of_usan_pairs.json
"""
import os
import json
import pandas as pd
from ast import literal_eval
from collections import Counter
from itertools import combinations


def match_drug_names(clinical_trials, drugs):
    """
    Matches drug names in clinical trial data to names in the drugs.csv file.
    Filters out common words using an ignore set.
    
    Args:
        clinical_trials (pd.DataFrame): Clinical trial data.
        drugs (pd.DataFrame): Drugs data from drugs.csv.
    
    Returns:
        list: List of dictionaries containing nct_id and matched drug names.
    """
    
    ignore_set = ('(seizure','prevention)', 'stimulating', 'hormone','(control)', '-', 'target', 
                  'controlled', 'infusion', 'system','(limited)','placebo', 'sugar', 'pill', '4%',
                  'gel', 'application', 'with', 'sham', 'microneedle', 'device', 'microneedle-facilitated',
                  'application', 'vein', 'soak', 'treated', 'with', 'lactated', 'ringers', 'solution', 
                  'hormone-related', 'protein', 'and', '2nd', 'phase', 'toxin', 'type-a', 'leaves', '5%',
                  'hydrochloride', 'continuation', 'omission', '1000mg', '2000mg', '250', 'mg', '500mg',
                  '50', 'iv', 'oral', 'antiplatelet', 'regimen', 'modification', 'or', 'methyl', 'recombinant',
                  'human', 'erythropoietin', 'cd', '2', 'administration', 'withdrawal', 'of', 'high-dose',
                  'regimens', 'use', '0.03%', 'topical', 'observational', 'antifungal', 'therapy', 'tissue',
                  'plasminogen', 'activator', 'normal', 'saline', 'adoptive', 'immunotherapy', 'single', 
                  'bolus', 'of', 'a', 'treatment', 'm,', 'reduced', 'nicotine', 'content', 'cigarettes', 'usual',
                  'patch', 'plus', 'gum/lozenge', 'acid', '100',  '200', 'cream', '0.3', '%', 'without', 'active',
                  'substance', '15%', 'inactivated', 'trivalent', 'influenza', 'vaccine', 'orally', 'everyday',
                  'hpv', 'vaccine', 'therapy.', '(identical', 'volume', 'of', 'normal', 'saline)', '0.9%', 'sodium',
                  'chloride', 'injectable', '1', 'gram', 'grams', 'target-controlled', 'anticoagulant', 'by', 
                  'physician', 'criteria', 'spray', 'nasal', '95%', 'pure',  'capsules', '200mg', 'tablet', 'group',
                  'adjuvant', 'perioperative', '(1-36)', 'fumarate', 'disoproxil', 'citrate', 's-1' ,'plus'
                 )

    # For every drug name identified in clinical trial data, match it to a drug name in drugs.csv
    for i, row in clinical_trials[
        clinical_trials.intervention_type == 'Drug'
    ].iterrows():
        lisst = []
        for j in [
            interv_name
            for interv_name in row['intervention_name'].split(' ')
            if interv_name not in ignore_set
        ]:
            try:
                lisst.append(
                    drugs[drugs.altLabel_list.str.contains(j)][
                        'itemLabel'
                    ].values[-1]
                )
            except:
                pass
        clinical_trials.loc[i, 'Drugs'] = str(lisst)

    output = clinical_trials[clinical_trials.intervention_type == 'Drug']
    output = output[output.Drugs != '[]']

    # JSON file of nct_id and drugs
    matched_drug_names = [
        {'nct_id': j.nct_id, 'drugs': literal_eval(j.Drugs)}
        for i, j in output[['nct_id', 'Drugs']].iterrows()
    ]

    with open(
        './output/matched_drug_names.json', 'w', encoding='utf-8'
    ) as file:
        json.dump(matched_drug_names, file)

    return matched_drug_names


def match_usan_codes(drug_names_file, usan_stems):
    """
     Associates drug names with their USAN prefixes/suffixes.

    Args:
        drug_names_file (str): Matched drug names from match_drug_names().
        usan_stems (pd.DataFrame): USAN stems data.

    Returns:
        list: List of dictionaries associating drugs with USAN codes and descriptions.
    """
    
    drug_names = pd.read_json(drug_names_file)

    # dealing with lists within lists
    drugs_list = list(
        set(drug for sublist in drug_names.drugs for drug in sublist)
    )

    # associate drug name with prefix and suffix
    drugs_usan = pd.DataFrame(drugs_list)
    drugs_usan_stem = pd.concat(
        [
            drugs_usan[0].str.contains(i)
            for i in usan_stems.stem.str.strip('-').dropna()
        ],
        axis=1,
    )
    drugs_usan_stem.columns = usan_stems.stem.str.strip('-').dropna().values
    drugs_usan = pd.DataFrame.join(drugs_usan, drugs_usan_stem)

    # index of drug names
    iter_index = drugs_usan[drugs_usan.iloc[:, 1:].sum(axis=1) != 0].index

    # determine USAN class
    def description(i):
        try:
            return {
                'description': usan_stems[
                    usan_stems.name.str.endswith(i) == True
                ]['definition'].values[0],
                'type': 'class',
            }
        except:
            try:
                return {
                    'description': usan_stems[
                        usan_stems.name.str.startswith(i) == True
                    ]['definition'].values[0],
                    'type': 'class',
                }
            except:
                return {
                    'description': usan_stems[
                        usan_stems.stem.str.endswith(i) == True
                    ]['definition'].values[0],
                    'type': 'subclass',
                }

    # JSON file identifying drug, usan_codes
    usan_json = [
        {
            'drug': drugs_usan.iloc[j, 0],
            'usan_codes': [
                description(i)
                for i in drugs_usan.iloc[j, 1:][drugs_usan.iloc[j, 1:]].index
            ],
        }
        for j in iter_index
    ]

    with open('./output/drugs_usan.json', 'w', encoding='utf-8') as file:
        json.dump(usan_json, file)

    return usan_json


def counts_of_trials_by_usan_class(drug_names_file, usan_json_file):
    """
    Associates clinical trials with USAN classes.

    Args:
        drug_names_file (str): matched drug names.
        usan_json_file (str): USAN data with associated drugs.

    Returns:
        list: List of USAN classes with associated trials.
    """
    
    drug_names = pd.read_json(drug_names_file).astype(str)
    usan_table = pd.read_json(usan_json_file)

    # associate trials to USAN class
    usan_table['trials'] = [
        list(
            drug_names[drug_names.drugs.str.contains(j.drug)]['nct_id'].values
        )
        for _, j in usan_table.iterrows()
    ]

    # add trails data into dict() of USAN class
    for i, j in usan_table.iterrows():
        for k in j.usan_codes:
            k['trials'] = usan_table['trials'][i]

    # JSON of USAN descriptions with associated trials info
    trials_by_usan = [
        k for _, j in usan_table.iterrows() for k in j.usan_codes
    ]

    with open('./output/trials_by_usan.json', 'w', encoding='utf-8') as file:
        json.dump(trials_by_usan, file)

    return trials_by_usan


def agg_counts_of_usan_pairs(trials_by_usan_file):
    """
    Aggregates counts of clinical trials associated for pairs of USAN classes.

    Args:
        trials_by_usan_file (str): USAN classes with associated trials.

    Returns:
        list: List of dictionaries containing trial counts for USAN class pairs.
    """

    trials = pd.read_json(trials_by_usan_file).astype(str)

    # ignore subclass
    trials = (
        trials[trials.type == 'class'].drop('type', axis=1).explode('trials')
    ).drop_duplicates()

    # list of unique paired combination of USAN descriptions for each trial
    pairs = trials.groupby(['trials']).agg(
        lambda g: list(set(combinations(sorted(g), 2)))
    )
    # pairs = pairs[pairs['description'].map(len) > 0]
    paired_count = pd.DataFrame(
        Counter(pairs.description.sum()).items(),
        columns=['USANclasspairs', 'tcount'],
    ).sort_values(
        'tcount', ascending=False # descending order
    )

    counts_of_usan_pairs = [
        {
            'description_1': j.USANclasspairs[0],
            'description_2': j.USANclasspairs[1],
            'trial_count': j.tcount,
        }
        for i, j in paired_count.iterrows()
    ]

    with open(
        './output/counts_of_usan_pairs.json', 'w', encoding='utf-8'
    ) as file:
        json.dump(counts_of_usan_pairs, file)

    return counts_of_usan_pairs


if __name__ == '__main__':
    # Import data
    clinical_trials = pd.read_json(
        'data/clinical_trials_2015.jsonl', lines=True
    )
    clinical_trials.intervention_name = (
        clinical_trials.intervention_name.str.lower()
    )

    drugs = pd.read_csv('data/drugs.csv').astype(str)
    drugs.altLabel_list = drugs.altLabel_list.str.lower()
    drugs.iloc[:, 1] = drugs.iloc[:, 1] + '|' + drugs.iloc[:, 0]

    usan_stems = pd.read_csv(
        'data/usan_stems.csv', on_bad_lines=print, engine='python'
    )

    # Results folder
    os.makedirs('./output', exist_ok=True)

    # Output files
    drug_names_file = './output/matched_drug_names.json'
    usan_json_file = './output/drugs_usan.json'
    trials_by_usan_file = './output/trials_by_usan.json'
    
    # Run tasks
    matched_drugs = match_drug_names(clinical_trials, drugs) # Q1
    usan_json = match_usan_codes(drug_names_file, usan_stems) #Q2
    trials_by_usan = counts_of_trials_by_usan_class(drug_names_file, usan_json_file) # Q3
    counts_of_usan_pairs = agg_counts_of_usan_pairs(trials_by_usan_file) # Q4
