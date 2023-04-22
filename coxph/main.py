import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

import matplotlib.pyplot as plt


def perform_coxph():
    # read in the clinical data
    clinical_df = pd.read_csv('data/lgg_clinical_with_rows.csv', index_col=0)
    clinical_df["index_for_merge"] = clinical_df.index

    # read in the mutation data
    mutation_df = pd.read_csv('data/lgg_mutation_with_rows.csv', index_col=0)
    mutation_df["index_for_merge"] = mutation_df["Tumor_Sample_Barcode"].apply(lambda x: x[:12].replace(
        "-", ".").lower())

    mutation_df = mutation_df[mutation_df['Hugo_Symbol'] == 'IDH1']

    # Merge clinical and mutation data based on sample IDs
    lgg_idh1_data = pd.merge(clinical_df, mutation_df,
                             left_on='index_for_merge', right_on='index_for_merge', how="outer")
    lgg_idh1_data.index = lgg_idh1_data["index_for_merge"]

    # silent mutation or not
    lgg_idh1_data['IDH1'] = lgg_idh1_data["Variant_Classification"].apply(
        lambda x: 1 if type(x) == str and x != "Silent" else 0)

    # specify gender: male = 1, female = 2
    lgg_idh1_data['gender'] = lgg_idh1_data["gender"].apply(
        lambda x: 1 if x == 'male' else 2)

    # specify race: white = 0, other = 1
    lgg_idh1_data['race'] = lgg_idh1_data["race"].apply(
        lambda x: 0 if x == 'white' else 1)

    # had radiation_therapy or not: therapy yes = 1, else =0
    lgg_idh1_data['radiation_therapy'] = lgg_idh1_data["radiation_therapy"].apply(
        lambda x: 1 if x == 'yes' else 0)

    # calculate duration_col
    lgg_idh1_data['duration_col'] = lgg_idh1_data["days_to_death"].add(
        lgg_idh1_data["days_to_last_followup"], fill_value=0)

    # filter patiants with 0 duration_col
    lgg_idh1_data = lgg_idh1_data[lgg_idh1_data['duration_col'] != 0]

    # Select relevant columns for analysis
    selected_columns = ['years_to_birth', 'vital_status', "gender", "radiation_therapy",
                        'duration_col', 'histological_type', 'race', 'IDH1']

    # 'karnofsky_performance_score' -- too many NaNs
    lgg_idh1_data = lgg_idh1_data[selected_columns]

    dummy_vars = pd.get_dummies(
        lgg_idh1_data[['histological_type']], drop_first=True)
    lgg_idh1_data = pd.concat([lgg_idh1_data, dummy_vars], axis=1)
    lgg_idh1_data.drop(['histological_type'], axis=1, inplace=True)

    # lgg_idh1_data.to_csv("data/processed.csv")
    # Drop rows with missing values
    lgg_idh1_data = lgg_idh1_data.dropna()
    lgg_idh1_data.to_csv("data/processed.csv")

    # Fit Cox proportional hazards model
    cph = CoxPHFitter()
    cph.fit(lgg_idh1_data, duration_col='duration_col',  # "days_to_death" or "days_to_last_followup"
            event_col='vital_status')

    # Print summary of fitted model
    # print(cph.summary)
    cph.summary.to_csv("result.csv")

    # PLOTING
    # plot the survival function
    # plot the survival function for the age variable
    fig, ax = plt.subplots()
    cph.plot_partial_effects_on_outcome(
        'years_to_birth', [-10, 0, 10, 20, 30, 40, 50, 60, 70, 100, 120, 150], cmap='coolwarm', ax=ax, plot_baseline=False)
    ax.set_xlabel('Overal survival (days)')
    ax.set_ylabel('Survival Probability')
    ax.set_title('Survival Function for Age')
    plt.tight_layout()

    # save the plot as a PNG file
    fig.savefig('survival_function_years.png')

    fig, ax = plt.subplots()
    cph.plot_partial_effects_on_outcome(
        'IDH1', [0, 1], cmap='coolwarm', ax=ax, plot_baseline=False)
    ax.set_xlabel('Overal survival (days)')
    ax.set_ylabel('Survival Probability')
    ax.set_title('Survival Function for IDH1')
    plt.tight_layout()

    # save the plot as a PNG file
    fig.savefig('survival_function_idh1.png')

    fig, ax = plt.subplots()
    cph.plot_partial_effects_on_outcome(
        'gender', [1, 2], cmap='coolwarm', ax=ax, plot_baseline=False)
    ax.set_xlabel('Overal survival (days)')
    ax.set_ylabel('Survival Probability')
    ax.set_title('Survival Function for gender')
    plt.tight_layout()

    # save the plot as a PNG file
    fig.savefig('survival_function_gender.png')


if __name__ == "__main__":
    perform_coxph()

    # try on 4-10 people
    # make a graph look normal

    # mutations -> each of them: how affects mortality. & together: how affects mortality
    # gene expression (?)

    # DNA variant
    # полиморфизм: >1% людей
    # мутация: <1% людей
    # HENG LI: читать (математика за ДНК, алгоримитка, etc)
    # microarray chip:
