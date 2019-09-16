import numpy as np
import pandas as pd
from scipy.stats import norm


class MixtureClassifier(object):
    def __init__(self, exp, mixture):
        """
        Classifies genes into normal and non-normal expression components

        :param pd.Series exp: Sample gene expression series in
                                     log2(TPM + 1)
        :param pd.DataFrame mixture: Mixture model parameters:
                                     Gene
                                     Normal_m - Normal mean
                                     Normal_std - Normal standard deviation
                                     Non_Normal_m - Non-Normal mean
                                     Non_Normal_std - Non-Normal standard deviation
        """
        # Only look at genes that are in the model
        idx = exp.index.intersection(mixture.index)
        self.mixture_df = mixture.loc[idx, :]
        self.exp_series = exp[idx]

    def fit(self):
        """
        Fits gene expression for each gene into one of two expression components.

        Assigns a 1 to a gene if the gene expression falls into the Non-Normal distribution and a
        0 otherwise. If either probability density returns NaN, assigns NaN to the gene in the
        classified series object.

        :return: classified genes
        :rtype: pd.Series
        """
        classified = pd.Series(index=self.mixture_df.index)
        for gene, row in self.mixture_df.iterrows():

            # Check if expression value is NaN
            if np.isnan(self.exp_series[gene]):
                raise ValueError('%s expression is NaN!' % gene)

            normal_norm = norm(row['Normal_m'], row['Normal_std'])
            not_normal_norm = norm(row['Non_Normal_m'], row['Non_Normal_std'])

            normal_lpdf = normal_norm.logpdf(self.exp_series[gene])
            not_normal_lpdf = not_normal_norm.logpdf(self.exp_series[gene])

            if not_normal_lpdf > normal_lpdf:
                classified[gene] = 1

            else:
                classified[gene] = 0

        return classified


if __name__ == '__main__':
    import sys

    if len(sys.argv) < 3:
        print("./classify-mixture-model.py exp_series mixture_df")
        exit(1)

    exp_df = pd.read_csv(sys.argv[1],
                         sep='\t',
                         index_col=0,
                         header=None)

    exp_series = pd.Series(index=exp_df.index,
                           data=exp_df[1])

    mixture_df = pd.read_csv(sys.argv[2],
                             sep='\t',
                             index_col=0)

    m = MixtureClassifier(exp_series, mixture_df)
    c = m.fit()
    print("Found %d non-normal genes" % sum(c))
