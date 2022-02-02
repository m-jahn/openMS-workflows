import pandas
import csv
import re
import numpy as np
import networkx as nx
from collections import defaultdict
from numpy import nansum, nanmean, array
from scipy import stats
from pyteomics import fasta


def fit_farms(probes, weight=0.5, mu=0, max_iter=1000, force_iter=False, min_noise=1e-4, fill_nan=1.0):
    '''
    Bayesian Factor Analysis for Proteomics Summarization
       A python translation of function "generateExprVal.method.farms" from Bioconductor FARMS.
       http://www.bioconductor.org/packages/release/bioc/html/farms.html

    Reference:
       Hochreiter S, Clevert D and Obermayer K (2006). A new summarization method for affymetrix probe level data.
       Bioinformatics, 22(8),
       http://bioinformatics.oxfordjournals.org/cgi/content/abstract/22/8/943.


    Inputs:
       probes: peptide abundance array (N peptides, M samples) in linear intensity scale
       weight: Hyperparameter value in the range of [0,1] which determines the influence of the prior.
       mu:     Hyperparameter value which allows to quantify different aspects of potential prior knowledge.
               A value near zero assumes that most genes do not contain a signal, and introduces a bias for
               loading matrix elements near zero.
    '''

    # normalize and transform X
    linear_readouts = pandas.DataFrame(probes)
    if fill_nan > 0:
        linear_readouts = linear_readouts.fillna(fill_nan)

    X = np.log2(linear_readouts).T
    X = X - np.nanmean(X, axis=0)
    xsd = np.nanstd(X, axis=0)
    xsd[xsd < min_noise] = 1.0
    X /= xsd
    X[~np.isfinite(X)] = 0

    n_samples, n_features = X.shape
    XX = np.cov(X.T, ddof=0)

    # positive definit
    XX = 0.5 * (XX + XX.T)
    XX[np.where(XX < 0)] = 0

    # robustness
    U,s,V = np.linalg.svd(XX)
    s[s < min_noise] = min_noise
    XX = U.dot(np.diag(s)).dot(V)

    # initiation
    L = np.sqrt(np.diag(XX) * 0.75)
    Psi = np.diag(XX) - L ** 2
    old_psi = Psi
    alpha = weight * n_features

    for i in range(max_iter):
        # E step
        PsiL = (1 / Psi) * L
        a = 1 + np.matrix(L) * np.matrix(PsiL).T
        bar = PsiL / a
        XXbeta = XX.dot(bar.T)
        Ezz = 1 - bar.dot(L) + bar.dot(XXbeta)

        # M step
        L = XXbeta.T / (Ezz + Psi * alpha)
        L = np.asarray(L)[0]
        Psi = np.diag(XX) - np.asarray(XXbeta)[0] * L + Psi * alpha * L * (mu - L)
        Psi[Psi < min_noise ** 2] = min_noise ** 2

        if abs(Psi - old_psi).max() / old_psi.max() < min_noise / 10:
            if not force_iter:
                break
        old_psi = Psi

    L *= np.sqrt(np.asarray(Ezz)[0]) # scale
    PsiL = (1 / Psi) * L
    SNR = 1. / (1 + np.matrix(L) * np.matrix(PsiL).T)

    weights = L / L.max()
    noise = SNR[0,0]
    return weights, noise


#weighted_average(loading, pep_abd, sampIx,
def weighted_average(weights, abundance, group_ix, w_cutoff=0.5, dof_loss=0, expMethod='geometric', impute_missing=0.001, nonmissing_threshold=0.01):
    assert expMethod in ('average', 'mean', 'geometric', 'gmean'), 'unknown averaging method!'
    expr = []
    aT = array(abundance)[weights > w_cutoff].T
    weights = weights[weights > w_cutoff]


    for ix in group_ix:
        if np.count_nonzero(np.nan_to_num(aT[ix])) < nonmissing_threshold * aT[ix].flatten().__len__():
            val = aT[ix]
            dof_loss +=  np.where(np.isnan(np.hstack(val)))[0].__len__()
            val[np.where(np.isnan(val))] = impute_missing
            aT[ix] = val
    if expMethod in ('average', 'mean'):
        pep_abd = np.multiply(aT, weights)
        for ix in group_ix:
            expr.append(np.nanmean(pep_abd[ix]))
    elif expMethod in ('geometric', 'gmean'):
        pep_abd = np.multiply(np.log2(aT), weights)
        for ix in group_ix:
            a_sum = 0
            w_sum = 0
            for a in pep_abd[ix]:
                a_sum += np.nansum(a)
                w_sum += np.sum(weights[np.isfinite(a)])
            expr.append(a_sum / w_sum)
        expr = 2 ** array(expr)
    pval = f_ANOVA(np.log2(aT.T), group_ix, np.log2(expr), dof_loss=dof_loss)
    return array(expr), pval


def sum_squares(peptide_abundances, sample_indeses, estimates=None):
    '''
    peptide_abundances:     n peptides, m samples
    sample_indeses:         k groups of sample indeses
    estimates:              k estimated abundances
    '''
    if estimates is None:
        est = np.array([abd[np.hstack(sample_indeses)] for abd in peptide_abundances])
        estimates = [np.nanmean(est) for _ in sample_indeses]

    residual = 0
    for i, est in zip(sample_indeses, estimates):
        res = peptide_abundances[:, i] - est
        residual += nansum(res ** 2)
    return residual


def f_ANOVA(peptide_abundances, sample_indeses, estimates, dof_loss=0):
    ss_total = sum_squares(peptide_abundances, sample_indeses)
    ss_resid = sum_squares(peptide_abundances, sample_indeses, estimates)
    dof1 = len(estimates) - 1
    dof2 = len([i for i in peptide_abundances.flatten() if np.isfinite(i)]) - len(estimates) - dof_loss
    f = ( (ss_total - ss_resid) / dof1) / (ss_resid / dof2)
    return stats.f.sf(f, dof1, dof2)



def read_peptide_protein_list(pep_file, pep_col=0, prot_col=None, header_rows=1, delimiter=','):
    '''
        param pep_file: text file with peptides and corresponding protein group
        outputs:
            a peptide set (I->L converted)
            an undirected graph (peptide <-> protein mapping)
    '''
    peps = set()
    g = nx.Graph()
    with open(pep_file) as fh:
        rd = csv.reader(fh, delimiter=delimiter)
        [rd.next() for _ in range(header_rows)] # skip header rows
        for row in rd:
            bseq = row[pep_col].replace('I', 'L').upper() # base amino acid chains with Ile -> Leu conversion
            peps.add(bseq)
            if not prot_col == None:
                g.add_edge(bseq, row[prot_col])
    return sorted(peps), g


def peptide_db_graph(peps, db, id_regex=None, dview = None):
    ''' search a set of peptides against a FASTA database  '''
    g = nx.Graph()
    prot_dict = dict()
    for header, seq, in fasta.read(db):
        seq = seq.replace('I', 'L').upper() # convert DB sequence I -> L
        prot_id = header.split()[0]
        if id_regex != None:
            find_id = re.findall(id_regex, header)
            if len(find_id) > 0:
                prot_id = find_id[0]
        prot_dict[prot_id] = seq

    def _map_seq(p):
        pairs = []
        for prot_id, seq, in prot_dict.items():
            if p in seq:
                pairs.append([p, prot_id])
        return pairs

    if dview:
        dview.push({'prot_dict': prot_dict})
        for ppps in dview.map_sync(_map_seq, list(peps)):
            if len(ppps):
                g.add_edges_from(ppps)
    else:
        for p in peps:
            ppps = _map_seq(p)
            if len(ppps):
                g.add_edges_from(ppps)
    return g


def parsimony_grouping(g, peps):
    '''
    Group peptides to proteins using the rule of parsimony
    g:  an undirected graph with peptide <-> protein as edges
    peps: the set of peptide sequences, nodes not listed in the peptide set are protein IDs.
    '''
    not_peps = set(g.nodes()) - set(peps)
    prot_groups = dict()
    for cc in nx.connected_component_subgraphs(g):
        in_group_peptides = set(cc.nodes()) - not_peps
        in_group_proteins = not_peps.intersection(cc.nodes())

        if len(in_group_proteins) == 1:
            prot_groups[in_group_proteins.pop()] = in_group_peptides
        elif len(in_group_proteins) > 1:
            reported = set()
            while len(in_group_proteins - reported) > 0:
                candidate_proteins = sorted(in_group_proteins - reported,
                                       key=lambda p:(len(set(cc[p].keys()) - reported), p),
                                       reverse=True)
                p = candidate_proteins[0]
                current_peps = set(cc[p].keys())
                plabel = [p]
                for i in range(1, len(candidate_proteins)):
                    _p = candidate_proteins[i]
                    _peps = set(cc[_p].keys())
                    if _peps == current_peps:
                        plabel.append(_p)
                    if len(_peps - current_peps) == 0:
                        reported.add(_p)

                plabel = ';'.join(sorted(plabel))
                if len(current_peps - reported) > 0:
                    prot_groups[plabel] = current_peps
                    reported = reported.union(current_peps)
                reported.add(p)
    return prot_groups


def protein_reader(df, sample_tag, group_tag):
    ''' Iterate peptide abundances from proteins in a pandas dataframe
        df: pandas dataframe containing peptide sequences, abundances and corresponding protein ID.
        sample_tag: list of sample IDs
        group_tag: column name of protein group
    '''
    df.sort(group_tag, inplace=True)
    values = df[samples].values
    tags = df[group_tag].values

    prot, peps = None, []
    for v, t in zip(values, tags):
        if prot == None:
            prot = t
        if t != prot:
            yield prot, array(peps)
            peps = []
            prot = t
        peps.append(v)
    if prot:
        yield prot, array(peps)


def read_peptide_matrix(fn, logScale=True, refSamples=''):
    df = pandas.read_csv(fn, index_col=0)
    df.index = [i.upper().replace('I', 'L') for i in df.index]
    return df


def protein_grouping(df, proteinDb):
    peptides = sorted(set(df.index))
    if not proteinDb:
        g = nx.Graph()
        for i, x in df.iterrows():
            for prot in x.values.astype('str')[0].split(';'):
                if len(prot) > 0:
                    g.add_edge(i, prot)
    else:
        g = peptide_db_graph(peptides, proteinDb)
    pg = parsimony_grouping(g, peptides)
    return pg


def zero_center_normalize(df, samples, logInput=False, method='median'):
    assert method in ('median', 'average', 'GMM'), 'Zero centering method has to be median, average or GMM!'
    if not logInput:
        # convert abundances to log2 scale
        df[samples] = df[samples].apply(np.log2)

    if method == 'average':
        norm_scale = np.nanmean(df[samples], axis=0)
    elif method == 'median':
        norm_scale = np.nanmedian(df[samples], axis=0)
    elif method == 'GMM':
        from sklearn.mixture import GMM
        # two-component Gaussian mixture model
        gmm = GMM(2)
        norm_scale = []
        for sp in samples:
            v = df[sp].values
            v = v[np.logical_not(np.isnan(v))]
            v = v[np.logical_not(np.isinf(v))]
            try:
                gmm.fit(np.matrix(v.values).T)
                vmean = gmm.means_[np.argmin(gmm.covars_)]
                norm_scale.append(vmean)
            except:
                norm_scale.append(np.nanmean(v))
        norm_scale = np.array(norm_scale)
    df[samples] = df[samples] - norm_scale
    return df



#=================================================
def main():
    import argparse, sys
    apars = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    apars.add_argument('-peptide',
                        default= 'peptides.csv',
                        help = 'Quantification result of peptides in CSV format.')
                        # The first column contains unique peptide sequences
                        # Missing values should be empty instead of zeros

    apars.add_argument('-db', nargs='?',
                        help='''Protein database in FASTA format.
                        If None, the peptide file must have protein ID(s) in the second column''')

    apars.add_argument('-samples', nargs='?',
                        help='''File of the sample list.
                        One run and its sample group per line, separated by tab.
                        If None, read from peptide file headings, then each run will be summarized as a group
                        ''')

    apars.add_argument('-log2', type=bool,
                        default=False,
                        help='Input abundances are in log scale (True) or linear scale (False)')

    apars.add_argument('-normalize', choices=['average', 'median', 'GMM', 'None'],
                        default='None', help='Method for sample-wise normalization.')
                        # Normalize input abundances (per sample) to zero-centered in log-scale
                        # Valid methods include: 'average', 'median' or 'GMM' (two-component Gaussian mixture model)
                        # If None (default), do not normalize

    # optional parameters
    apars.add_argument('-mu', type=float,
                        default=0.1, help='Hyperparameter mu')
                        # Hyperparameter mu of the FARMS algorithm: prior knowledge of the expected loading.

    apars.add_argument('-alpha', type=float,
                        default=0.1, help='Hyperparameter weight of prior probability')
                        # Hyperparameter weight of the FARMS algorithm: weight of prior probability in EM calculation.

    apars.add_argument('-ref',
                        default='average', help='Names of reference sample groups (separated by semicolon)')
                        # If average (default), calculate average peptide abundance as the reference.
                        # Otherwise, keep abundance values as is.

    apars.add_argument('-min_sample', type=int,
                        default=1, help='Minimum number of samples peptides needed to be quantified in')
                        # Peptides quantified in less than the minimum number will be discarded


    apars.add_argument('-pep_unique', type=bool,
                        default=False, help='Use unique peptides only')

    apars.add_argument('-nonmissing', type=float,
                        default=0.01,
                        help= '''Minimum fraction of non-missing values in the group.
                                Imputing missing values if non-missing fraction is lower than the threshold.'
                                ''')

    apars.add_argument('-cutoff_weight', type=float,
                        default=0.5, help='Peptides weighted lower than the cutoff will be excluded')

    apars.add_argument('-fast', type=bool,
                        default=True, help='Allow early termination in EM calculation.')
                        # When the noise is sufficiently small (i.e. force 1000 times iterations)

    apars.add_argument('-out', type=argparse.FileType('w'),
                        default=sys.stdout, help='Path to output file (writing in TSV format).')

    apars.add_argument('-out_type', default='protein', help='Output file rows: Either \'protein\' or \'peptide\'')

    args = apars.parse_args()


    #------------------------------------------------
    df = read_peptide_matrix(args.peptide, args.log2) #log2 'Input abundances are in log scale (True) or linear scale (False)' Default: False
    print("Abundance matrix loaded: %d peptides" % len(df.index))

    #samples_f='/home/vitalv/cyano_dataset_20170714/diffacto/samples_f_original_ctrls.csv'
    if not args.samples:
        samples = [s for s in df.columns]
        if args.db == None:
            samples.pop(0)
            print(args.db)
        groups = samples
    else:
        samples = []
        groups = []
        with open(args.samples) as fh:
            for line in fh.readlines():
                try:
                    _s, _g = line.rstrip().split('\t')
                    samples.append(_s)
                    groups.append(_g)
                except ValueError:
                    pass

    # per sample normalization of peptide abundances # Normalize input abundances (per sample) to zero-centered in log-scale. Deafult: None
    if not args.normalize == 'None':
        df = zero_center_normalize(df, samples, logInput=args.log2, method=args.normalize)
        args.log2 = True

    # select reference runs if specified
    ref_samples = []
    if args.ref:
        for r in args.ref.split(';'):
            for i in range(len(groups)):
                if groups[i] == r:
                    ref_samples.append(i)
    ref_samples = [samples[i] for i in ref_samples]

    print("Number of runs: %d" % len(samples))
    print("Reference runs: %s" % len(ref_samples))
    print(*ref_samples, sep='\t')

    # sample grouping
    group_names = [i for i in sorted(set(groups), key=lambda k: "{0:0>50}".format(k)) if i not in args.ref.split(';')]
    sampIx = np.array([[j for j in range(len(groups)) if groups[j] == i] for i in group_names])
    print("Number of sample groups: %d" % len(group_names))

    # coverage filtering # Peptides quantified in less than the minimum number (min_sample) will be discarded.
    # min_sample Minimum number of samples peptides needed to be quantified in. Default: 1
    df = df[[np.count_nonzero(np.nan_to_num(v)) >= int(args.min_sample) for v in df[samples].values]]


    # protein grouping
    pg = protein_grouping(df, args.db) #db='/home/vitalv/database/Synechocystis_PCC6803_protein_sequences.fasta'
    print("Number of protein groups: %d" % len(pg.keys()))
    # print(*sorted(pg.keys()), sep='\n')

    # reversed mapping (peptide to protein group) for checking peptide uniqueness.
    pep2prot = defaultdict(list)
    for prot_ids, bseqs in pg.items():
        for s in bseqs:
            pep2prot[s] += prot_ids.split()

    # use unique peptides #Deafualt False
    if args.pep_unique:
        df = df[[len(pep2prot[p]) == 1 for p in df.index]]


    #------------------------------------------------------------------------------
    # perform differential analysis
    if not args.out_type or args.out_type == 'protein':
        print(*['Protein', 'N.Pept', 'Q.Pept', 'S/N', 'Pvalue'] + group_names, sep='\t', file=args.out)
    elif args.out_type == 'peptide':
        print(*['peptide'] + ['protein'] + ['weight'] + group_names , sep='\t', file=args.out)

    for prot in sorted(pg.keys()):
        peps = pg[prot]
        dx = df.ix[[p for p in peps if p in df.index]]
        pep_count = len(dx)
        dof_loss = 0

        if args.log2: #Input abundances are in log scale (True) or linear scale (False)' Default: False
            dx[samples] = 2 ** dx[samples]

        pep_abd = dx[samples].values


        if len(ref_samples):
            # use reference for peptide abundance re-scaling
            reference_abundance = dx[ref_samples].mean(axis=1).fillna(np.nanmean(dx[samples])).values
            pep_abd = (pep_abd.T / reference_abundance).T

        elif args.ref.lower() == 'average':
            # use average abundances for re-scaling
            reference_abundance = dx[samples].mean(axis=1).values
            pep_abd = (pep_abd.T / reference_abundance).T
        else:
            # do not re-scale
            pass


        if pep_count == 1: # single peptide group
            loading = array([1 for _ in dx.index])
            noise = 1.0
            #continue # not report
        elif pep_count > 1:
            loading, noise = fit_farms(pep_abd,mu=args.mu,weight=args.alpha,max_iter=1000,force_iter=not args.fast,fill_nan=1.0)
        else:
            continue


        sn = 10 * np.log10((1-noise)/noise)
        qc_passed = len(loading[loading>args.cutoff_weight]) #'Peptides weighted lower than the cutoff will be excluded . Default 0.5
        if args.ref == 'average':
            dof_loss = qc_passed

        protein_summary_group, pval = weighted_average(loading, pep_abd, sampIx,
                                        w_cutoff=args.cutoff_weight,
                                        dof_loss=dof_loss,
                                        expMethod='gmean',
                                        impute_missing=np.nanmin(pep_abd)/2,
                                        nonmissing_threshold=args.nonmissing)

        if not args.out_type or args.out_type == 'protein':
            print(prot, pep_count, qc_passed, sn, pval, *np.log2(protein_summary_group), sep='\t', file=args.out)
        elif args.out_type == 'peptide':
            l=0
            for i, row in dx.iterrows(): #uncomment for all peptides output
                print(i, prot, loading[l], *row.values,  sep='\t', file=args.out)
                l+=1




if __name__ == '__main__':
    main()
