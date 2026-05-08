import re, tempfile
from dataclasses import dataclass
from pathlib import Path
import pandas as pd
from biorazer.design.basic import Entry, Library


@dataclass
class SingleBindcraft(Entry):

    pass


@dataclass
class BatchBindcraft(Library):

    entry_type = SingleBindcraft

    @staticmethod
    def explode_df(df: pd.DataFrame):
        df = df.reset_index(drop=True)
        sample_is = set()
        shared_cols = []
        col_suffixes = set()
        for col in df.columns:
            matched = re.match(r"(\d+)_(.*)", col)
            if matched:
                sample_is.add(int(matched.group(1)))
                col_suffixes.add(matched.group(2))
            else:
                if not col.startswith("Average"):
                    shared_cols.append(col)
        sample_is = list(sample_is)
        sample_is.sort()
        col_suffixes = list(col_suffixes)
        all_data = []
        for row_idx in df.index:
            for sample in sample_is:
                row_data = dict()
                row_data.update(df.loc[row_idx, shared_cols].to_dict())
                row_data["sample"] = sample
                for s in col_suffixes:
                    row_data[s] = df.at[row_idx, f"{sample}_{s}"]
                all_data.append(row_data)
        samples = pd.DataFrame(all_data)
        selector = pd.notna(samples["pLDDT"])
        samples = samples[selector].reset_index(drop=True)
        return samples
