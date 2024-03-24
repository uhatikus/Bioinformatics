import os
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter

import matplotlib.pyplot as plt


class CoxPH():
    def __init__(self, clinical_file_path: str, mutation_file_path: str, results_path: str, skip_mutations: bool = False) -> None:
        # read in the clinical data
        self.clinical_df = pd.read_csv(
            clinical_file_path, index_col=9, sep=";")
        # self.clinical_df["index_for_merge"] = self.clinical_df.index

        self.skip_mutations = skip_mutations

        # read in the mutation data
        if not self.skip_mutations:
            self.mutation_df = pd.read_csv(
                mutation_file_path, index_col=0)
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

        if not self.skip_mutations:
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
                lambda x: 1 if x == 'MALE' else 2)

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
                            'duration_col'] + covariats

        preprocessed_clinical_data = preprocessed_clinical_data[selected_columns]

        if 'histological_type' in covariats:
            dummy_vars = pd.get_dummies(
                # histological_type
                preprocessed_clinical_data[['histological_type']])
            preprocessed_clinical_data = pd.concat(
                [preprocessed_clinical_data, dummy_vars], axis=1)
            preprocessed_clinical_data.drop(
                ['histological_type'], axis=1, inplace=True)
            preprocessed_clinical_data.drop(
                ["histological_type_oligoastrocytoma"], axis=1, inplace=True)

        return preprocessed_clinical_data

    def preprocess_mutation_data_with_gene(self, gene: str):
        cur_mutation_df = self.mutation_df[self.mutation_df['Hugo_Symbol'] == gene]
        cur_mutation_df = cur_mutation_df[[
            "Variant_Classification", "index_for_merge"]]

        cur_mutation_df.dropna(inplace=True)
        # silent mutation or not
        # Variant_Classification:
        # 1. Missence
        # 2. PTV = Frame_shift, Nonsence, Splice site
        # Others -- not used
        cur_mutation_df[gene] = cur_mutation_df["Variant_Classification"].apply(
            lambda x: 1 if x != "Silent" else 0)
        cur_mutation_df.drop(
            columns=["Variant_Classification"], inplace=True)

        return cur_mutation_df

    def plot_survival_function(self, data: pd.DataFrame, target_column: str, options: list[(str, any)], cutoff: str = ""):
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

        fig.savefig(self.results_path +
                    f'survival_function_{target_column}_{cutoff}.png')

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

    def analize_genes_together(self, genes: list[str], covariats: list[str] = None):
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
            if gene in cur_data.columns:
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

    def analize_genes_together_with_ATRX(self, genes: list[str], covariats: list[str] = None, target_gene: str = "ATRX"):
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
            if gene in cur_data.columns:
                cur_data[gene].fillna(0, inplace=True)

            # drop nans
            cur_data.dropna(inplace=True)

        print("ready for CoxPH")

        for cov in cur_data.columns:
            if cov == target_gene or cov == 'duration_col' or cov == 'vital_status':
                continue
            cur_data_i = cur_data.drop(columns=[cov])

            cph = CoxPHFitter()
            cph.fit(cur_data_i, duration_col='duration_col',
                    event_col='vital_status')
            # Save summary of fitted model
            cph.summary.to_csv(self.results_path +
                               f"result_all_genes_without_{cov}.csv")

            if target_gene in cph.summary[cph.summary['p'] < 0.05].index:
                print(f"If {cov} in covariates, then {target_gene} has p < 0.05")
            else:
                print(
                    f"If {cov} in covariates, then {target_gene} DOES NOT HAVE p < 0.05")

    def valid_covariats(self, covariats: list[str]):
        if covariats is None:
            return self.default_covariats
        else:
            for cov in covariats:
                if cov not in self.default_covariats:
                    raise Exception(
                        f"{cov} is not possible covariant for the current dataset")
        return covariats

    def analize_genes(self, genes: list[str], covariats: list[str] = None, target_gene: str = None):
        covariats = self.valid_covariats(covariats)

        good_genes = []
        bad_genes = []
        for gene in genes:
            # len(genes) == 1
            try:
                self.analize_gene(gene, covariats, plot_flag=False)
                good_genes.append(gene)
            except:
                # print(f"Gene {gene} gives an error => cannot be calculated")
                bad_genes.append(gene)
        print(f"Unprocessed gense: {bad_genes}")

        self.analize_genes_together(good_genes, covariats)

        if target_gene is not None:
            self.analize_genes_together_with_ATRX(
                good_genes, covariats, target_gene)

    def analyze_genes_sum_effect(self, genes: list[str], covariats: list[str] = None):
        covariats = self.valid_covariats(covariats)
        cur_data = self.preprocess_clinical_data_using_covariats(
            covariats)
        present_genes = []
        bad_genes = []
        for gene in genes:
            clean_mutation_data_for_current_gene = self.preprocess_mutation_data_with_gene(
                gene)
            if clean_mutation_data_for_current_gene.empty:
                # print(
                #     f'gene {gene} does not present in the current coxph dataset')
                bad_genes.append(gene)
                continue
            # Merge clinical and mutation data for current gene
            cur_data = pd.merge(cur_data, clean_mutation_data_for_current_gene,
                                left_on='index_for_merge', right_on='index_for_merge', how="outer")

            cur_data.index = cur_data["index_for_merge"]
            cur_data.drop(columns=["index_for_merge"], inplace=True)
            if gene in cur_data.columns:
                cur_data[gene].fillna(0, inplace=True)
                present_genes.append(gene)
            else:
                print(
                    f'gene {gene} does not present in the current coxph dataset !!!!!!')

        cur_data["mutation effect"] = cur_data[present_genes].sum(axis=1)
        cur_data.dropna(inplace=True)
        cur_data.to_csv(self.results_path+"data_sum_of_mutations.csv")
        cur_data.drop(columns=present_genes, inplace=True)

        print("ready for CoxPH")

        cph = CoxPHFitter()
        cph.fit(cur_data, duration_col='duration_col',
                event_col='vital_status')

        print("CoxPH fitted")
        # Save summary of fitted model
        cph.summary.to_csv(self.results_path+f"result_sum_of_mutations.csv")

        self.plot_survival_function(cur_data, target_column="mutation effect", options=[
            (2, f"2 mutations"), (1, f"1 mutation"), (0, f"No mutation")])
        # print(f"bad genes: {bad_genes}")

    def analyze_by_two_groups(self, defined_group: list[str], covariats: list[str] = None):
        covariats = self.valid_covariats(covariats)
        cur_data = self.preprocess_clinical_data_using_covariats(
            covariats)
        defined_group = [x[:12].replace("-", ".").lower()
                         for x in defined_group]
        cur_data['group_1'] = cur_data.index.isin(defined_group).astype(int)

        cur_data.dropna(inplace=True)
        cur_data.drop(columns=["index_for_merge"], inplace=True)

        print(cur_data)
        cph = CoxPHFitter()
        cph.fit(cur_data, duration_col='duration_col',
                event_col='vital_status')

        cph.summary.to_csv(self.results_path+f"group_1.csv")

        self.plot_survival_function(cur_data, target_column="group_1", options=[
            (1, f"in group"), (0, f"Not in group")])
        # print(f"bad genes: {bad_genes}")

    def analyze_by_two_defined_groups(self, group1: list[str], group2: list[str], cutoff: str, covariats: list[str] = None):
        covariats = self.valid_covariats(covariats)
        cur_data = self.preprocess_clinical_data_using_covariats(
            covariats)

        cur_data['group'] = (1) * cur_data.index.isin(group1).astype(
            int) + (-1) * cur_data.index.isin(group2).astype(int)

        # cur_data['group2'] = cur_data.index.isin(group2).astype(int)

        cur_data.dropna(inplace=True)
        # cur_data.drop(columns=["index_for_merge"], inplace=True)

        # print(cur_data['duration_col'])
        # print(cur_data['vital_status'])
        # print(cur_data['group'])
        # print(cur_data.columns)
        cph = CoxPHFitter()
        cph.fit(cur_data, duration_col='duration_col',
                event_col='vital_status')

        cph.summary.to_csv(self.results_path +
                           f"matched_vs_unmatched_"+cutoff+".csv")

        self.plot_survival_function(cur_data, target_column="group", options=[
            (-1, f"unmatched"), (1, f"matched"), (0, f"others")], cutoff=cutoff)
        # print(f"bad genes: {bad_genes}")


if __name__ == "__main__":
    coxph = CoxPH(clinical_file_path='coxph/data/lgg/lgg_clinical_with_rows.csv',
                  mutation_file_path='coxph/data/lgg/lgg_mutation_with_rows.csv',
                  results_path='coxph/results/some_results/')
