# set snakemake bin 
snakemake := snakemake

demo:
	-mkdir -p $@

.PHONY: install run tscc clean debug

install: enviroment.yaml
	conda env create --file $<

debug : demo Snakefile
	cp config.yaml $</config.yaml
	cp pre.py $</pre.py
	cp $(word 2,$^) $</Snakefile
	cd $< && ${snakemake} --use-conda \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete --skip-script-cleanup


run : demo Snakefile
	cp config.yaml $</config.yaml
	cp pre.py $</pre.py
	cp $(word 2,$^) $</Snakefile
	cd $< && ${snakemake} \
     -c 1 -p --snakefile Snakefile -R --rerun-incomplete

tscc: demo Snakefile
	cp config.yaml $</config.yaml
	cp $(word 2,$^) $</Snakefile
	cp -R pbs_profile/. $</pbs_profile
	cd $< && ${snakemake} --config content=bye\
     -c 1 -p --snakefile Snakefile -R \
     --latency-wait 20 \
     --rerun-incomplete --profile pbs_profile

clean:
	-rm -rf demo
