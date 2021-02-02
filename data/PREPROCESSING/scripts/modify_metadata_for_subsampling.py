from argparse import ArgumentParser
from collections import defaultdict
from collections import ChainMap
import pandas as pd
from Bio import SeqIO
from random import choice

def get_lineages(f):
    return [x for x in pd.read_csv(f)['lineage'].unique() if x!='None']

def set_bel_lineage(df,lineages):
    df['be_lineage'] = df.apply(lambda x: int(any([x['country']=='Belgium',
                                        x['pangolin_lineage'] in lineages])), axis=1)
    df['be_lineage'] = df['be_lineage'].apply(lambda x: "yes" if x==1 else 'no')
    
def case_tiers(f):
    casedf = pd.read_csv(f).query('continent == "Europe"')
    casedf['country'] = casedf['country'].str.replace("Czechia",'Czech Republic')
    # remove rows tracking EUR and UK+EUR cases
    countries = [x for x in casedf['country'].unique() if "(total)" not in x]
    # count total number of cases by country
    counts = (casedf.groupby("country")['weekly_count']
                .sum().sort_values(ascending=False)
                .reset_index()
                .query('country in @countries')
                .set_index('country')['weekly_count'])
    total_cases = counts.sum()
    cum_perc = counts.cumsum()/total_cases
    case_tier = pd.DataFrame()
    case_tier['tot_cases'] = counts
    case_tier['perc'] = cum_perc
    return case_tier

def eur_case_tier(row):
    if row['region']!='Europe':
        return "NA"
    elif row['country'] in tier1:
        return 'T1'
    elif row['country'] in tier2:
        return 'T2'
    else:
        return 'T3'

def set_eur_case_tier(df):
    df['case_tier'] = df.apply(eur_case_tier, axis=1)
    
def load_seqs(f,meta):
    seqs = {}
    for r in SeqIO.parse(f,'fasta'):
        #partially trim utr
        seqs[r.id] = str(r.seq[100:-100])
    seqnames = list(seqs.keys())
    seq_meta = meta.query('strain in @seqnames').copy()
    seq_meta['month'] = seq_meta['date'].apply(lambda x:x.split("-")[1])
    return seqs, seq_meta

def dedup(group, seqdict):
    genome_hash = defaultdict(list)
    for s in group:
        genome = seqdict[s]
        genome_hash[genome].append(s)
    keep = [choice(cluster) for cluster in genome_hash.values()]        
    return {s:('yes' if s in keep else 'no') for s in group},list(genome_hash.values())

def dedup_identical(s,seqdict):
    dedup_clusters =  s.apply(lambda x:dedup(x,seqdict))
    all_clusts = dedup_clusters.apply(lambda x:x[-1][0]).to_frame()
    dedup_dict = [x[0] for x in dedup_clusters.values]
    return ChainMap(*dedup_dict),all_clusts
    

def date_window(dt):
    try:
        dt = pd.to_datetime(dt)
        if dt <= pd.to_datetime("2020-03-31"):
            return '2020Q1'
        elif dt in pd.date_range('2020-04-01','2020-06-30'):
            return '2020Q2'
        elif dt in pd.date_range('2020-07-01','2020-09-30'):
            return '2020Q3'
        elif dt in pd.date_range("2020-10-01","2020-12-31"):
            return '2020Q4'
        else:
            assert dt >= pd.to_datetime("2021-01-01")
            return '2021'
    except:
        return "NA"
    
if __name__=="__main__":

    # Argument parser for the command line
    parser = ArgumentParser("Modifies default metadatafile for custom subsampling")
    parser.add_argument("-m", "--metadata", type=str, required=True, help="Name of the nextstrain metadata file to modify.")
    parser.add_argument("-p", "--pangolin", type=str, required=True, help="Name of the pangolin clade assignment output CSV.")
    parser.add_argument("-c", "--cases", type=str, required=True, help="Name of the ECDC case counts csv file")
    parser.add_argument("-s", "--seqs", type=str, required=True, help="Fasta alignment of all Belgian sequences")
    parser.add_argument("-o", "--out", type=str,required=True, help="output file name")
    args = parser.parse_args()
    
    nextmeta = pd.read_csv(args.metadata, low_memory=False,sep='\t')
    print('... Adding lineage filter ....')
    set_bel_lineage(nextmeta, get_lineages(args.pangolin))
    cases = case_tiers(args.cases)
    tier1 = cases.loc[cases['perc']<0.5].index.values
    tier2 = tier2 = cases.query('perc >=0.5').query('perc <0.9').index.values
    print('... Adding case tier for EUR countries ...')
    set_eur_case_tier(nextmeta)
    seqdict, seq_meta = load_seqs(args.seqs, nextmeta)
    grouped_seqs = seq_meta.groupby(["division","month"])['strain'].apply(lambda x:list(x))
    print('... Deduplicating identical BE sequences...')
    dedup_dict,dedup_clusters = dedup_identical(grouped_seqs, seqdict)
    dedup_clusters.to_csv("BE_identical_seq_clusters.csv")
    nextmeta['dedup'] = nextmeta['strain'].map(dedup_dict).fillna("NA")
    nextmeta['time_window'] = nextmeta['date'].apply(date_window)

    nextmeta.to_csv(args.out,sep='\t',index=False)
    
    
    
        
    
                  
    
    
    