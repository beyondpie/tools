# set snakemake bin 
snakemake := snakemake

demo:
	-mkdir -p $@

.PHONY: install run tscc clean

install: enviroment.yaml
	conda env create --file $<

run : demo Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && ${snakemake} \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete

tscc: demo Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R profile/. $</profile
	cd $< && ${snakemake} --config content=bye\
     -c 1 -p --snakefile Snakefile -R \
     --latency-wait 20 \
     --rerun-incomplete --profile profile

clean:
	-rm -rf demo
