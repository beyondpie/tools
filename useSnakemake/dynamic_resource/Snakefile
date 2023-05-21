## optional config file
configfile: "config.yaml"

content = config["content"]
samples = config["samples"]
## all the output files you want to have
rule all:
    input:
        expand("flag/pre_{s}.done", s = samples),
        expand("flag/first_{s}.done", s = samples)

## then set up the rules about how to generate them
rule pre:
    output:
        # s will be infered based on rule all output
        # here it will be a or b.
        # snakemake will run s=a and s=b in parallel if possible
        # touch will be automatically  generate the flag file
        # once the rule is done.
        touch("flag/pre_{s}.done")
    log:
        "log/{s}.log"
    shell:
        """
        # wildcards.s to get the a / b
        echo "pre:" {content} {wildcards.s} 2> {log}
        """

# snakemake will know that it depends on the output of pre
# then the rule will run after pre
rule first:
    input:
        "flag/pre_{s}.done"
    output:
        touch("flag/first_{s}.done")
    log:
        "log/{s}.log"
    shell:
        """
        echo "first:" {content} {wildcards.s} 2> {log}
        """
