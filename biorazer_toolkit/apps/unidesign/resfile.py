import pandas as pd

ALL_AA = "ACDEFGHIKLMNPQRSTVWY"
AA_EXCEPT_CYS = "ADEFGHIKLMNPQRSTVWY"


class RESFILE:

    def __init__(self, sites_design_df: pd.DataFrame = None):
        """

        Parameters
        ----------
        sites_design_df : pd.DataFrame
            DataFrame with columns: ['chain', 'res_id', 'res_name1"]
        """
        if sites_design_df is None:
            self.sites_design_df = pd.DataFrame(
                columns=["chain", "res_id", "res_name1"]
            )
        else:
            required_columns = {"chain", "res_id", "res_name1"}
            if not required_columns.issubset(sites_design_df.columns):
                raise ValueError(
                    f"sites_design_df must contain columns: {required_columns}"
                )
            self.sites_design_df = sites_design_df

    def set_design(self, chain_id, res_id, res_name1):
        """
        Set design for a specific residue.

        Parameters
        ----------
        chain_id : str
            Chain identifier
        res_id : int
            Residue number
        res_name1 : str
            String consisting of one-letter code of the amino acid to design to
        """
        mask = (self.sites_design_df["chain"] == chain_id) & (
            self.sites_design_df["res_id"] == res_id
        )
        if not mask.any():
            new_row = pd.DataFrame(
                {"chain": [chain_id], "res_id": [res_id], "res_name1": [res_name1]}
            )
            self.sites_design_df = pd.concat(
                [self.sites_design_df, new_row], ignore_index=True
            )
        else:
            self.sites_design_df.loc[mask, "res_name1"] = res_name1

    def write(self, file):
        """

        Parameters
        ----------
        file : str
            Path to output resfile (a .txt file)
        """
        with open(file, "w") as f:
            f.write("SITES_DESIGN_START\n")
            for _, row in self.sites_design_df.iterrows():
                chain = row["chain"]
                residue_number = row["res_id"]
                amino_acid = row["res_name1"]
                f.write(f"{chain} {residue_number:>6d} {amino_acid}\n")
            f.write("SITES_DESIGN_END\n")
