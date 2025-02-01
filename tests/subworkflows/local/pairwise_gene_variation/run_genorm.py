import pandas as pd
import numpy as np
import sys

file = sys.argv[1]
# Expression data for three control genes.
counts = pd.read_parquet(file)
counts.set_index("ensembl_gene_id", inplace=True)
counts = counts.T.replace(0, 1e-8)


def _m_numpy(gene_expression: np.ndarray) -> np.ndarray:
    """Internal control gene-stability measure `M`.

    Computes Eq. (4) in Ref. [1].

    [1]: Vandesompele, Jo, et al. "Accurate normalization of real-time quantitative
    RT-PCR data by geometric averaging of multiple internal control genes." Genome
    biology 3.7 (2002): 1-12.
    """

    if not (gene_expression > 0).all():
        raise ValueError(
            "Expression domain error: not all expression data are strictly positive!"
        )

    a = gene_expression
    # Eq. (2): A_{jk}^{(i)} = log_2 (a_{ij} / a_{ik})
    A = np.log2(np.einsum("ij,ik->ijk", a, 1 / a))
    # Eq. (3)
    V = np.std(A, axis=0)
    # Eq. (4) N.B., Since V_{j=k} is zero, we can simply ignore it since it does not
    # contribute to calculation.
    n = V.shape[1]
    return np.sum(V, axis=1) / (n - 1)


def m_measure(gene_expression):
    m_values = _m_numpy(gene_expression.to_numpy())
    return pd.Series(m_values, index=gene_expression.columns)


print(m_measure(counts).sort_values())
