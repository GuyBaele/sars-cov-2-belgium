import pandas as pd
from argparse import ArgumentParser

def read_metadata(f):
	print(f"reading metadata file {f}...")
	df = pd.read_csv(f,sep='\t',low_memory=False)
	df['fulldate'] = df['date'].apply(lambda x: "XX" not in x)
	df = df.query("fulldate == True").copy()
	df['date'] = pd.to_datetime(df['date'])
	df['week'] = df['date'].apply(lambda x:f"{x.year}-{x.week}")
	return df

def subsample(df, n, lineage, lineage_prop,max_BE_lin):
	if lineage is not None:
		# sample selected lineage
		lin_seqs = df.query(f"pango_lineage == '{lineage}'")
		# check if number of sequences from lineage satisfies lineage_prop asked
		if len(lin_seqs)/n <= lineage_prop:
			print(f"number of {lineage} sequences available can't satisfy proportion asked")
			print(f"subsample size: {n}")
			print(f"{lineage} sequences available: {len(lin_seqs)}")
			print(f"lineage proportion asked: {lineage_prop}")
			print(f"Subsample size will be adjusted to keep lineage proportion")
			n = len(lin_seqs)/lineage_prop

		num_lineage_seqs = int(n*lineage_prop)
		# include all lineage sequences available if rescaling was necessary
		print(f"Subsampling {lineage} sequences...")
		if num_lineage_seqs == len(lin_seqs):
			LIN_SAMPLE = lin_seqs
		else:
			num_BE_sample = min([max_BE_lin,len(lin_seqs.query('country=="Belgium"'))])
			num_nonBE_sample = num_lineage_seqs - num_BE_sample
			print(f"  Subsampling Belgian {lineage} sequences...")
			BE_VOC_SAMPLE = generic_sampling(lin_seqs.query('country=="Belgium"'),num_BE_sample,"division")
			print(f"  Subsampling non-Belgian {lineage} sequences...")
			WORLD_VOC_SAMPLE = generic_sampling(lin_seqs.query('country!="Belgium"'),num_nonBE_sample,"country")
			LIN_SAMPLE = pd.concat([BE_VOC_SAMPLE,WORLD_VOC_SAMPLE])
		print(f"Subsampling non-{lineage} sequences...")
		NON_LIN_SAMPLE = generic_sampling(df.query(f"pango_lineage != '{lineage}'"),n-num_lineage_seqs,"country")
		ALL_SAMPLED = pd.concat([LIN_SAMPLE, NON_LIN_SAMPLE])
		bel_voc_selected = len(LIN_SAMPLE.query('country=="Belgium"'))
		print(f"{len(ALL_SAMPLED)} sequences sampled from a total of {len(df)}")
		print(f" * {len(LIN_SAMPLE)} sequences of lineage {lineage}")
		print(f" * {bel_voc_selected} Belgian sequences of lineage {lineage}")

		return ALL_SAMPLED

	else:
		print("No lineage to focus sampling was selected")

	return None

def generic_sampling(df,num,geo_level):
	# start by sampling at least one from each geo_level:
	initial_sample = (df.groupby(geo_level)
					  .apply(lambda x: x.sample(1))
					  .reset_index(drop=True))
	initial_ids = initial_sample['gisaid_epi_isl'].values
	if num <= len(initial_ids):
		return initial_sample
	else:
		missing_seqs = int(num - len(initial_ids))
		remaining_seqs = df.query('gisaid_epi_isl not in @initial_ids')
		# sample one seqs per location and week
		sample2 = (remaining_seqs
				   .groupby(["division",'week'])
				   .apply(lambda x:x.sample(1))
				   .reset_index(drop=True))
		if len(sample2) >= missing_seqs:
			sample2 = sample2.sample(missing_seqs)
			return pd.concat([initial_sample,sample2])
		else:
			sample2_ids = sample2['gisaid_epi_isl'].values
			missing2 = num - len(initial_ids) - len(sample2_ids)
			remaining2 = remaining_seqs.query("gisaid_epi_isl not in @sample2_ids")
			if missing2 <= len(remaining2):
				# sample3 = remaining2.sample(missing2)
				num_locs = len(remaining2[geo_level].unique())
				sample3_1 = (remaining2.groupby(geo_level)
							 .apply(lambda x: x.sample(int(missing2/num_locs),replace=True))
							 .reset_index(drop=True)
							 .drop_duplicates())
				ids_3_1 = sample3_1['gisaid_epi_isl'].values
				sample3_2 = (remaining2.query("gisaid_epi_isl not in @ids_3_1")
						   .sample(missing2 - len(ids_3_1))
						   .reset_index(drop=True))
				sample3 = pd.concat([sample3_1,sample3_2])
			else:
				sample3 = remaining2
			return pd.concat([initial_sample,sample2,sample3])

if __name__ == "__main__":
	parser = ArgumentParser()
	parser.add_argument("--metadata",required=True,help='BEAST XML file')
	parser.add_argument("--n",type=int,default=2500,help='output file name')
	parser.add_argument("--lineage",default=None,type=str,help='lineage to focus on')
	parser.add_argument("--lineage_prop",default=None,type=float,help='proportion of sequences to sample lineage selected')
	parser.add_argument("--max_BE_lin",default=300,type=float,help='maximum number of Belgian sequences from lineage selected')
	parser.add_argument("--out",required=True)
	args = parser.parse_args()
	if args.lineage is not None:
		args.lineage_prop = 0.5

	print(f"Focus sampling on {args.lineage}")
	df = read_metadata(args.metadata)
	sampled_seqs = subsample(df, args.n, args.lineage, args.lineage_prop, args.max_BE_lin)
	if sampled_seqs is not None:
		sampled_seqs.to_csv(args.out,index=False)
