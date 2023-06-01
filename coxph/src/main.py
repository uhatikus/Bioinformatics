from coxph import CoxPH

root_folder = 'coxph/'


def test_one_gene(gene):
    coxph = CoxPH(clinical_file_path=f'{root_folder}data/lgg/lgg_clinical_with_rows.csv',
                  mutation_file_path=f'{root_folder}data/lgg/lgg_mutation_with_rows.csv',
                  results_path=f'results/one_gene/{gene}/')
    coxph.analize_genes([gene])


def test_unique_genes_from_dataset():
    MIN_NUMBER_OF_ENTRIES = 20
    coxph = CoxPH(clinical_file_path=f'{root_folder}data/lgg/lgg_clinical_with_rows.csv',
                  mutation_file_path=f'{root_folder}data/lgg/lgg_mutation_with_rows.csv',
                  results_path=f'{root_folder}results/unique_genes_results/')

    gene_counts = coxph.mutation_df["Hugo_Symbol"].value_counts()
    mask = coxph.mutation_df['Hugo_Symbol'].isin(
        gene_counts[gene_counts < MIN_NUMBER_OF_ENTRIES].index)
    unique_genes = set(
        coxph.mutation_df[~mask]["Hugo_Symbol"].unique())
    coxph.analize_genes(unique_genes)


def test_unique_genes_from_dataset_and_analyze_gene(target_gene="ATRX"):
    MIN_NUMBER_OF_ENTRIES = 20
    coxph = CoxPH(clinical_file_path=f'{root_folder}data/lgg/lgg_clinical_with_rows.csv',
                  mutation_file_path=f'{root_folder}data/lgg/lgg_mutation_with_rows.csv',
                  results_path=f'{root_folder}results/unique_genes_with_{target_gene}_results/')

    gene_counts = coxph.mutation_df["Hugo_Symbol"].value_counts()
    mask = coxph.mutation_df['Hugo_Symbol'].isin(
        gene_counts[gene_counts < MIN_NUMBER_OF_ENTRIES].index)
    unique_genes = set(
        coxph.mutation_df[~mask]["Hugo_Symbol"].unique())

    coxph.analize_genes(unique_genes, target_gene=target_gene)


def test_genes_from_file(file_path: str):  # 'lgg_genes_module.txt'
    coxph = CoxPH(clinical_file_path=f'{root_folder}data/lgg/lgg_clinical_with_rows.csv',
                  mutation_file_path=f'{root_folder}data/lgg/lgg_mutation_with_rows.csv',
                  results_path=f'{root_folder}results/{file_path.split("/")[-1].split(".")[0]}_results/')

    with open(file_path, 'r') as file:
        genes = set([gene.strip() for gene in file.readlines()])
        coxph.analize_genes(genes)


def display_candidates(genes, genes_from_file, unique_genes_from_datatset):
    print()
    print("Candidates:")
    for g_i in (genes_from_file-genes):
        f = True
        for g_j in unique_genes_from_datatset:
            if g_i in g_j and g_i != g_j:
                if f:
                    print(f"\n{g_i}: {g_j}", end="")
                    f = False
                else:
                    print(f", {g_j}", end="")
    print()


def test_sum_mutations():
    coxph = CoxPH(clinical_file_path=f'{root_folder}data/lgg/lgg_clinical_with_rows.csv',
                  mutation_file_path=f'{root_folder}data/lgg/lgg_mutation_with_rows.csv',
                  results_path=f'{root_folder}results/sum_of_mutations_results/')

    unique_genes_from_datatset = set(
        coxph.mutation_df["Hugo_Symbol"].unique())
    with open(f'{root_folder}given_files/lgg_genes_module.txt', 'r') as file:
        genes_from_file = set([gene.strip() for gene in file.readlines()])
    print(
        f"Not present in the current dataset: {genes_from_file-unique_genes_from_datatset}")
    genes = unique_genes_from_datatset & genes_from_file
    print(f"Present in the current dataset: {genes}")
    coxph.analyze_genes_sum_effect(list(genes))

    display_candidates(genes, genes_from_file, unique_genes_from_datatset)


def skcm_analysis():
    coxph = CoxPH(clinical_file_path=f'{root_folder}data/skcm/skcm_clinical_with_rows.csv',
                  mutation_file_path=f'{root_folder}data/skcm/skcm_mutation_with_rows.csv',
                  results_path=f'{root_folder}results/skcm_results/', skip_mutations=True)

    with open(f'{root_folder}given_files/skcm_atm_mis_ptv_indiv.txt', 'r') as file:
        people = [gene.strip() for gene in file.readlines()]

    coxph.analyze_by_two_groups(defined_group=people, covariats=[
        "gender", 'years_to_birth'])  # , 'stage', 'grade'


if __name__ == "__main__":
    print("let's start")
    skcm_analysis()


# updated TODO:
# 1. ATRX and NF1 literature study.

# TODO:
# ATRX "research": drop one or another covariate and see the reuslt (p-value for ATRX): https://www.frontiersin.org/articles/10.3389/fonc.2017.00236/full


# If years_to_birth in covariates, then ATRX DOES NOT HAVE p < 0.05

# If gender in covariates, then ATRX has p < 0.05
# If radiation_therapy in covariates, then ATRX has p < 0.05
# If race in covariates, then ATRX has p < 0.05
# If histological_type_astrocytoma in covariates, then ATRX has p < 0.05
# If histological_type_oligodendroglioma in covariates, then ATRX has p < 0.05
# If MUC16 in covariates, then ATRX has p < 0.05
# If NOTCH1 in covariates, then ATRX has p < 0.05
# If PIK3CA in covariates, then ATRX has p < 0.05

# If NF1 in covariates, then ATRX DOES NOT HAVE p < 0.05

# If FUBP1 in covariates, then ATRX has p < 0.05
# If IDH1 in covariates, then ATRX has p < 0.05
# If TP53 in covariates, then ATRX has p < 0.05
# If TTN in covariates, then ATRX has p < 0.05
# If CIC in covariates, then ATRX has p < 0.05

# TODO:
# Unite all the mutations from the list. Comparison: No mutations from the list, some mutations from the list. More mutation => more effect.
# MYC: MYCBP, MYCBP2
# GEM: GEMIN6, GEMIN5
# APP: PAPPA2, APPL1, PAPPA, APPBP2, TRAPPC8, TRAPPC11
# AQP1: AQP10
# KRT19: KRT19P2
# GPR3: GPR31, GPR39, GPR32, GPR35

# pattern of gene expression (find similar) -> consider only missence mutations
# e-mail with github (uhatikus)

# TO CONSIDER
# 1. Mutations (done)
# 2. Gene expression (done)
# 3. Tumor micro environment (missence mutations burden)  - на сколько ткань отличима от врождённой (???)
# https://gnomad.broadinstitute.org/


# TODO
# wait for gene expression
# wait for code
# wait for people id (positive = mutation exists) -> compare people (coxph), do diff expression

# TODO:
# TCGA on python -- search - does exist?:
# - https://www.google.com/search?q=TCGA+on+python&oq=TCGA+on+python&aqs=chrome..69i57.262j0j7&sourceid=chrome&ie=UTF-8
# - https://github.com/arahuja/pytcga
# - https://github.com/hammerlab/pygdc

# something that exists in python, but not in R (like hail: https://github.com/hail-is/hail)


# ______________
# theory of 3 (somatic) mutations
# if 1 "born" mutation => 2 somatic mutations (cancer earlier)
# gene ATM.
# is there a signature of gene expression if there is a mutation in gene ATM?


# TODO:
# https://pubmed.ncbi.nlm.nih.gov/35740680/

# TODO:
# 1. T-test
# 2. log_2 (average(vector group 1))/average(vector group 2))
# some genes will have low p-value
# 3. pathway: http://www.pantherdb.org/

# микроокружение опохули влияет или нет?

# https://www.nature.com/articles/s41467-022-31436-8
