import xml.etree.ElementTree as ET

def root_xpath_text(r, p):
    tmp = r.find(p)
    assert tmp is not None, f"XPath failed: {p}"
    return tmp.text

def root_xpath_attribute_value(r, p, a):
    tmp = r.find(p)
    assert tmp is not None, f"XPath failed: {p}"
    assert a in tmp.attrib, f"{a} missing."
    return tmp.attrib[a]

tree = ET.parse("config.xml")
root = tree.getroot()

time_series_clean_csv = root_xpath_text(root, ".//results/intermediate/timeSeries")
data_plot_png = root_xpath_text(root, ".//results/figures/dataPlot")
data_plot_rds = data_plot_png.replace("png", "rds")
disaster_txt = root_xpath_text(root, ".//results/intermediate/disasterStrings")
fasta = root_xpath_text(root, ".//results/intermediate/sequences")
present_rds = root_xpath_text(root, ".//results/intermediate/present")
beast_bin = "lib/beast/bin/beast"
combined_fig_png = root_xpath_text(root, ".//combinedEverything")

# The following pulls out the MCMC XML file and then constructs the filepath of the logger file that this will
mcmc_xml = root_xpath_text(root, ".//results/intermediate/beastXML")
mcmc_tree = ET.parse(mcmc_xml)
mcmc_root = mcmc_tree.getroot()
maybe_mcmc_log = root_xpath_attribute_value(mcmc_root, ".//logger[@id=\"tracelog\"]", "fileName")
if "$(filebase)" in maybe_mcmc_log:
    mcmc_log = maybe_mcmc_log.replace("$(filebase)", mcmc_xml.split("/")[-1].replace(".xml", ""))
else:
    mcmc_log = maybe_mcmc_log


# This is the main entry point for the snakemake pipeline.
rule all:
    input:
        beast_bin,
        time_series_clean_csv,
        data_plot_png,
        disaster_txt,
        fasta,
        present_rds,
        mcmc_log,
        combined_fig_png


rule plot_combined_results:
    input:
        "R/postprocessing-1.R",
        "config.xml",
        data_plot_rds,
        present_rds,
        mcmc_xml,
        mcmc_log
    output:
        combined_fig_png
    shell:
        """
        Rscript {input[0]}
        """


rule run_mcmc:
    input:
        mcmc_xml,
        beast_bin
    output:
        mcmc_log
    shell:
        """
        ./lib/beast/bin/beast -seed 1 -overwrite {input[0]}
        """


rule download_beast_and_tracer:
    input:
        "download-beast.sh"
    output:
        beast_bin
    shell:
        "bash download-beast.sh"


rule make_time_series_and_data_plot:
    input:
        "R/preprocessing-1.R",
        "config.xml"
    output:
        time_series_clean_csv,
        data_plot_png,
        data_plot_rds
    shell:
        """
        Rscript {input[0]}
        """


rule make_fasta_and_disaster_files:
    input:
        "R/preprocessing-2.R",
        "config.xml",
        time_series_clean_csv
    output:
        disaster_txt,
        fasta,
        present_rds
    shell:
        """
        Rscript {input[0]}
        """
