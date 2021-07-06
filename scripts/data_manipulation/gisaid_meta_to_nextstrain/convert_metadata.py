import pandas as pd
import argparse
from Bio import SeqIO

def get_seq_len(s):
    seqdict={}
    for r in SeqIO.parse(s,'fasta'):
        seqdict[r.description] = len(r.seq)
    return seqdict

def fix_seq_names(s):
    n = f"updated_{s}"
    records = []
    for r in SeqIO.parse(s,'fasta'):
        r.id = sf(str(r.id.split("|")[0].split(" ")[0]))
        r.description = ""
        records.append(r)
    SeqIO.write(records,n,"fasta")

def sf(s):
    if s.startswith("hCoV-19"):
        s = s[8:]
    return s


def fill_incomplete_date(dt):
    assert type(dt)==str,'date must be a string'
    resolution = len(dt.split('-'))
    if resolution == 3:
        return dt
    else:
        return dt + '-XX'*(3-resolution)

def format_columns(df,seq_lens):
    df=df.copy()
    nextstrain_cols = ['strain','virus','gisaid_epi_isl','genbank_accession','date',
        'region','country', 'division', 'location', 'region_exposure', 'country_exposure',
        'division_exposure', 'segment', 'length', 'host', 'age', 'sex', 'Nextstrain_clade',
        'pangolin_lineage','GISAID_clade', 'originating_lab','submitting_lab', 'authors',
        'url', 'title','paper_url', 'date_submitted','purpose_of_sequencing']

    # df['strain'] = df[['Virus name','Accession ID','Collection date']].apply(lambda x:"|".join(x),axis=1)
    df['strain'] = df['Virus name'].apply(sf)
    df['virus'] = 'ncov'
    df['gisaid_epi_isl'] = df['Accession ID']
    df['genbank_accession'] = '?'
    df['date'] = df['Collection date'].apply(fill_incomplete_date)
    df['region'] = df['Location'].apply(lambda x: x.split(' / ')[0] if x.count(' / ') >=1 else '?')
    df['country'] = df['Location'].apply(lambda x: x.split(' / ')[1] if x.count(' / ') >=1 else '?')
    df['division'] = df['Location'].apply(lambda x: x.split(' / ')[2] if x.count(' / ') >=2 else '?')
    df['location'] = df['Location'].apply(lambda x: x.split(' / ')[3] if x.count(' / ') >=3 else '?')
    df['region_exposure'] = '?'
    df['country_exposure'] = '?'
    df['division_exposure'] = '?'
    df['segment'] = 'genome'
    df['length'] = df['strain'].map(seq_lens)
    df['host'] = df['Host']
    df['age'] = df['Patient age']
    df['sex'] = df['Gender']
    df['Nextstrain_clade'] = '?'
    df['pangolin_lineage'] = df['Lineage']
    df['GISAID_clade'] = df['Clade']
    df['originating_lab'] = '?'
    df['submitting_lab'] = '?'
    df['authors'] = '?'
    df['url'] = 'https://www.gisaid.org'
    df['title'] = '?'
    df['paper_url'] = '?'
    df['date_submitted'] = '?'
    df['purpose_of_sequencing'] = '?'
    return df[nextstrain_cols].copy()

def read_gisaid_meta(seqMeta,patientMeta):
    df1 = pd.read_csv(seqMeta,sep='\t')
    df2 = pd.read_csv(patientMeta,sep='\t')
    shared_cols = list(set(df1.columns.to_list()).intersection( set(df2.columns.to_list())))
    return pd.merge(left=df1,right=df2,on=shared_cols)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True, help="FASTA file of sequences")
    parser.add_argument("--seqMeta", required=True, type=str, help="GISAID sequence metadata")
    parser.add_argument("--patientMeta", required=True, help="GISAID patient metadata")
    parser.add_argument("--nextMeta", required=True, help="Nextmeta tsv")
    parser.add_argument("--outFile",required=True,help='output file')

    args = parser.parse_args()

    gisaid_df = read_gisaid_meta(args.seqMeta,args.patientMeta)
    seq_lens = get_seq_len(args.fasta)
    nextMeta = pd.read_csv(args.nextMeta,sep='\t')
    new_nextmeta = format_columns(gisaid_df,seq_lens)
    updated_meta = pd.concat([nextMeta,new_nextmeta],verify_integrity=True,ignore_index=True)
    updated_meta.to_csv(args.outFile,sep='\t',index=False)
    fix_seq_names(args.fasta)
