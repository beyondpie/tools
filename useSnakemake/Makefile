demo:
	-mkdir -p $@

.PHONY: run tscc

run : demo Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete

tscc: test_dir Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && snakemake --config content=bye\
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --profile profile
