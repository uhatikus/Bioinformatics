
import pandas as pd
from lifelines import CoxPHFitter

if __name__ == "__main__":

    # read in the clinical data
    clinical_df = pd.read_csv('data/lgg_clinical_with_rows.csv', index_col=0)
    clinical_df["index_for_merge"] = clinical_df.index

    # read in the mutation data
    mutation_df = pd.read_csv('data/lgg_mutation_with_rows.csv', index_col=0)
    mutation_df["index_for_merge"] = mutation_df["Tumor_Sample_Barcode"].apply(lambda x: x[:12].replace(
        "-", ".").lower())
    # Merge clinical and mutation data based on sample IDs
    lgg_idh1_data = pd.merge(clinical_df, mutation_df,
                             left_on='index_for_merge', right_on='index_for_merge')
    lgg_idh1_data.index = lgg_idh1_data["index_for_merge"]

    # filter IDH1 Hugo_Symbol
    lgg_idh1_data = lgg_idh1_data[lgg_idh1_data['Hugo_Symbol'] == 'IDH1']

    # silent mutation or not
    # lgg_idh1_data['harmful'] = lgg_idh1_data["Variant_Classification"].apply(
    #     lambda x: 1 if x != 'Silent' else 0)

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
                        'duration_col', 'histological_type', 'race']  # 'harmful', "Hugo_Symbol"

    # 'karnofsky_performance_score' -- too many NaNs
    lgg_idh1_data = lgg_idh1_data[selected_columns]

    print(list(pd.factorize(
        lgg_idh1_data['histological_type'])[1]))
    # 'astrocytoma' : 0, 'oligoastrocytoma' : 1, 'oligodendroglioma' : 2
    lgg_idh1_data['histological_type'] = pd.factorize(
        lgg_idh1_data['histological_type'])[0]

    # Drop rows with missing values
    lgg_idh1_data = lgg_idh1_data.dropna()
    lgg_idh1_data.to_csv("data/processed.csv")

    # Fit Cox proportional hazards model
    cph = CoxPHFitter()
    cph.fit(lgg_idh1_data, duration_col='duration_col',  # "days_to_death" or "days_to_last_followup"
            event_col='vital_status')

    # Print summary of fitted model
    print(cph.summary)
    cph.summary.to_csv("result.csv")
