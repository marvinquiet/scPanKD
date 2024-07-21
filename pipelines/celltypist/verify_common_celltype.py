import fire
import pandas as pd


def print_celltypes(
    data_path: str = 'datasets_common_celltype/ProjecTIL_EGAS00001004809_CD8T_trimmed_metadata.csv',
    col: str = 'curated_celltype',
):
    data = pd.read_csv(data_path)
    print(f'data from {data_path}')
    print(data[col].value_counts())


if __name__ == '__main__':
    fire.Fire()