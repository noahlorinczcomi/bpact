import yaml
import argparse
import pymc as pm
import numpy as np
import pandas as pd


def load_config(path: str) -> dict:
    with open(path, "r") as f:
        config = yaml.safe_load(f)
    return config


def main(args: argparse.ArgumentNamespace) -> None:
    """
    1. Load GWAS
    2. Clean GWAS cols
    3. Load gene-gene correlation matrix (on the fly)
    4. 
    """
    cfg = args.config

    return None


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="BPACT CLI: Estimate number of causal genes using GWAS sumstats"
    )
    p.add_argument(
        "-i", "--config", type=str, required=True, help="Full path to config ``.yaml`` file"
    )
    args = p.parse_args()
    main(args)
