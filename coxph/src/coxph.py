import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter

import matplotlib.pyplot as plt


class CoxPH():
    def __init__(self, clinical_file_path: str, mutation_file_path: str, results_path: str) -> None:
        # read in the clinical data
        self.clinical_df = pd.read_csv(
            'data/lgg_clinical_with_rows.csv', index_col=0)
        self.clinical_df["index_for_merge"] = self.clinical_df.index

        # read in the mutation data
        self.mutation_df = pd.read_csv(
            'data/lgg_mutation_with_rows.csv', index_col=0)
        self.mutation_df["index_for_merge"] = self.mutation_df["Tumor_Sample_Barcode"].apply(lambda x: x[:12].replace(
            "-", ".").lower())

        self.validate_data()

        self.results_path = results_path
        if not os.path.exists(results_path):
            os.mkdir(results_path)

        self.default_covariats = ['years_to_birth', "gender", "radiation_therapy",
                                  'histological_type', 'race']
        self.update_default_covariats()

    def validate_data(self):
        def validate(df: pd.DataFrame, required_colums: list[str]):
            for col in required_colums:
                if col not in df.columns:
                    raise Exception(
                        f"clinical data does not have {col}. Please provide clinical data with {col}")

        validate(self.clinical_df, ['vital_status',
                                    'days_to_death', 'days_to_last_followup'])

        validate(self.mutation_df, ['Variant_Classification'])

    def update_default_covariats(self):
        new_default_covariats = []
        for cov in self.default_covariats:
            if cov in self.clinical_df.columns:
                new_default_covariats.append(cov)
        self.default_covariats = new_default_covariats

    def preprocess_clinical_data_using_covariats(self, covariats: list[str]) -> pd.DataFrame:
        preprocessed_clinical_data = self.clinical_df.copy()

        # specify gender: male = 1, female = 2
        if 'gender' in covariats:
            preprocessed_clinical_data['gender'] = preprocessed_clinical_data["gender"].apply(
                lambda x: 1 if x == 'male' else 2)

        # specify race: white = 0, other = 1
        if 'race' in covariats:
            preprocessed_clinical_data['race'] = preprocessed_clinical_data["race"].apply(
                lambda x: 0 if x == 'white' else 1)

        # had radiation_therapy or not: therapy yes = 1, else =0
        if 'radiation_therapy' in covariats:
            preprocessed_clinical_data['radiation_therapy'] = preprocessed_clinical_data["radiation_therapy"].apply(
                lambda x: 1 if x == 'yes' else 0)

        # calculate duration_col
        preprocessed_clinical_data['duration_col'] = preprocessed_clinical_data["days_to_death"].add(
            preprocessed_clinical_data["days_to_last_followup"], fill_value=0)

        # filter patiants with 0 duration_col
        preprocessed_clinical_data = preprocessed_clinical_data[
            preprocessed_clinical_data['duration_col'] != 0]

        # Select relevant columns for analysis
        selected_columns = ['vital_status',
                            'duration_col', 'index_for_merge'] + covariats

        preprocessed_clinical_data = preprocessed_clinical_data[selected_columns]

        if 'histological_type' in covariats:
            dummy_vars = pd.get_dummies(
                # histological_type
                preprocessed_clinical_data[['histological_type']], drop_first=True)
            preprocessed_clinical_data = pd.concat(
                [preprocessed_clinical_data, dummy_vars], axis=1)
            preprocessed_clinical_data.drop(
                ['histological_type'], axis=1, inplace=True)

        return preprocessed_clinical_data

    def preprocess_mutation_data_with_gene(self, gene: str):
        cur_mutation_df = self.mutation_df[self.mutation_df['Hugo_Symbol'] == gene]
        cur_mutation_df = cur_mutation_df[[
            "Variant_Classification", "index_for_merge"]]

        cur_mutation_df.dropna(inplace=True)
        # silent mutation or not
        cur_mutation_df[gene] = cur_mutation_df["Variant_Classification"].apply(
            lambda x: 1 if x != "Silent" else 0)
        cur_mutation_df.drop(
            columns=["Variant_Classification"], inplace=True)

        return cur_mutation_df

    def plot_survival_function(self, data: pd.DataFrame, target_column: str, options: list[(str, any)]):
        kmf = KaplanMeierFitter()

        T = data["duration_col"]
        E = data["vital_status"]
        fig, ax = plt.subplots()

        sum_m = 0

        for i, option in enumerate(options):
            m = (data[target_column] == option[0])
            kmf.fit(durations=T[m], event_observed=E[m],
                    label=option[1])
            # at_risk_counts gives some statistics (for the first option only)
            kmf.plot_survival_function(ax=ax, at_risk_counts=(i == 0))
            sum_m = m + sum_m

        if len(sum_m.unique()) != 1 or sum_m.unique() != 1:
            print(
                f"Not options are considered. Possible options:{data[target_column].unique()}")

        ax.set_xlabel('Overal survival (days)')
        ax.set_ylabel('Survival Probability')
        plt.title(f'Survival Function for {target_column}')

        fig.savefig(self.results_path+f'survival_function_{target_column}.png')

    def analize_gene(self, gene: str, covariats: list[str], plot_flag: bool = False):
        clean_clinical_data = self.preprocess_clinical_data_using_covariats(
            covariats)
        clean_mutation_data_for_current_gene = self.preprocess_mutation_data_with_gene(
            gene)

        # Merge clinical and mutation data for current gene
        cur_data = pd.merge(clean_clinical_data, clean_mutation_data_for_current_gene,
                            left_on='index_for_merge', right_on='index_for_merge', how="outer")
        cur_data.index = cur_data["index_for_merge"]
        cur_data.drop(columns=["index_for_merge"], inplace=True)

        # of no gene mutation for the current gene and patient, then 0
        cur_data[gene].fillna(0, inplace=True)

        # drop nans
        cur_data.dropna(inplace=True)

        cph = CoxPHFitter()
        cph.fit(cur_data, duration_col='duration_col',
                event_col='vital_status')

        # Save summary of fitted model
        cph.summary.to_csv(self.results_path+f"result_{gene}.csv")

        if plot_flag:
            self.plot_survival_function(cur_data, target_column=gene, options=[
                (1, f"{gene}"), (0, f"No {gene}")])
            self.plot_survival_function(cur_data, target_column="gender", options=[
                (1, "male"), (2, "female")])
            self.plot_survival_function(cur_data, target_column="race", options=[
                (0, "white"), (1, "others")])
            self.plot_survival_function(cur_data, target_column="radiation_therapy", options=[
                (1, "yes"), (0, "no"), ])

    def analize_genes_togather(self, genes: list[str], covariats: list[str] = None):
        cur_data = self.preprocess_clinical_data_using_covariats(
            covariats)
        for gene in genes:
            clean_mutation_data_for_current_gene = self.preprocess_mutation_data_with_gene(
                gene)

            # Merge clinical and mutation data for current gene
            cur_data = pd.merge(cur_data, clean_mutation_data_for_current_gene,
                                left_on='index_for_merge', right_on='index_for_merge', how="outer")
            cur_data.index = cur_data["index_for_merge"]
            cur_data.drop(columns=["index_for_merge"], inplace=True)

            # of no gene mutation for the current gene and patient, then 0
            cur_data[gene].fillna(0, inplace=True)

            # drop nans
            cur_data.dropna(inplace=True)

        print("ready for CoxPH")

        cph = CoxPHFitter()
        cph.fit(cur_data, duration_col='duration_col',
                event_col='vital_status')

        print("CoxPH fitted")
        # Save summary of fitted model
        cph.summary.to_csv(self.results_path+f"result_all_genes.csv")

        cph.summary[cph.summary['p'] < 0.05].to_csv(
            self.results_path+f"significant_covariants_and_genes.csv")

    def analize_genes(self, genes: list[str], covariats: list[str] = None):
        if covariats is None:
            covariats = self.default_covariats
        else:
            for cov in covariats:
                if cov not in self.default_covariats:
                    raise Exception(
                        f"{cov} is not possible covariant for the current dataset")
        for gene in genes:
            # len(genes) == 1
            self.analize_gene(gene, covariats, plot_flag=True)

        self.analize_genes_togather(genes, covariats)


if __name__ == "__main__":
    coxph = CoxPH(clinical_file_path='../data/lgg_clinical_with_rows.csv',
                  mutation_file_path='../data/lgg_mutation_with_rows.csv',
                  results_path='results/')

    gene_counts = coxph.mutation_df["Hugo_Symbol"].value_counts()
    mask = coxph.mutation_df['Hugo_Symbol'].isin(
        gene_counts[gene_counts < 20].index)
    unique_genes = coxph.mutation_df[~mask]["Hugo_Symbol"].unique()
    genes = ["IDH1", "NOS1"]
    genes = ["DNAH8"]
    coxph.analize_genes(unique_genes)
