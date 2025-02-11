import xml.etree.ElementTree as ET
import os

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

# --------------------------------------------------------------------
# The following pulls out the MCMC XML file and then constructs the
# filepath of the logger file that this will

mcmc_xml = root_xpath_text(root, ".//results/intermediate/beastXML")
mcmc_tree = ET.parse(mcmc_xml)
mcmc_root = mcmc_tree.getroot()
maybe_mcmc_log = root_xpath_attribute_value(mcmc_root, ".//logger[@id=\"tracelog\"]", "fileName")
if "$(filebase)" in maybe_mcmc_log:
    mcmc_log = maybe_mcmc_log.replace("$(filebase)", mcmc_xml.split("/")[-1].replace(".xml", ""))
else:
    mcmc_log = maybe_mcmc_log
# --------------------------------------------------------------------
# The following pulls out the data that is needed for the subsampling
# experiment.

ss_disaster_json = root_xpath_text(root, ".//subsampling/intermediate/disasterStringsJSON")
ss_disaster_png = root_xpath_text(root, ".//subsampling/figures/disasterPlot")
ss_xml_dir = root_xpath_text(root, ".//subsampling/intermediate/newXMLDir")
# TODO The XML filepaths should really be constructed from the
# configuration data, but because that is a bit of a pain, I'm just
# going to hardcode it for now.
ss_orig_xml = "out/subsampling-experiment/xml/timtam-2023-10-15-original.xml"
ss_orig_log = "out/timtam-2023-10-15-original.log"
ss_mcmc_xml_a = "out/subsampling-experiment/xml/timtam-2023-10-15-subsample-0p66.xml"
ss_mcmc_log_a = "out/timtam-2023-10-15-subsample-0p66.log"
ss_mcmc_xml_b = "out/subsampling-experiment/xml/timtam-2023-10-15-subsample-0p33.xml"
ss_mcmc_log_b = "out/timtam-2023-10-15-subsample-0p33.log"
# --------------------------------------------------------------------


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
        combined_fig_png,
        # Subsampling Experiment
        ss_disaster_json,
        ss_disaster_png,
        ss_disaster_png.replace("png", "svg"),
        ss_disaster_png.replace("png", "rds"),
        ss_orig_xml,
        ss_mcmc_xml_a,
        ss_mcmc_xml_b,
        ss_orig_log,
        ss_mcmc_log_a,
        ss_mcmc_log_b,
        "out/subsampling-experiment/summary-plot-r0.png",
        "out/subsampling-experiment/summary-plot-r0.svg",
        "out/subsampling-experiment/summary-plot-r0.rds",
        "out/subsampling-experiment/summary-plot-historysizes.png",
        "out/subsampling-experiment/summary-plot-historysizes.svg",
        "out/subsampling-experiment/summary-plot-historysizes.rds",
        "out/subsampling-experiment/summary-plot-obs_props.png",
        "out/subsampling-experiment/summary-plot-obs_props.svg",
        "out/subsampling-experiment/summary-plot-obs_props.rds",
        "out/manuscript/subsampling-experiment-combined-r0-timeseries.png"


rule plot_manuscript_subsampling_combined:
    input:
        "R/postprocessing-6-subsampling.R",
        "out/subsampling-experiment/summary-plot-r0.rds",
        "out/subsampling-experiment/summary-plot-historysizes.rds",
        ss_disaster_png.replace("png", "rds"),
    output:
        "out/manuscript/subsampling-experiment-combined-r0-timeseries.png"
    shell:
        """
        Rscript {input[0]}
        """


rule plot_subsampling_obs_props_comparison:
    input:
        ss_orig_log,
        ss_mcmc_log_a,
        ss_mcmc_log_b,
        rscript="R/postprocessing-5-subsampling.R",
    output:
        "out/subsampling-experiment/summary-plot-obs_props.png",
        "out/subsampling-experiment/summary-plot-obs_props.svg",
        "out/subsampling-experiment/summary-plot-obs_props.rds"
    shell:
        """
        Rscript {input.rscript}
        """


rule plot_subsampling_historysize_comparison:
    input:
        ss_orig_log,
        ss_mcmc_log_a,
        ss_mcmc_log_b,
        rscript="R/postprocessing-4-subsampling.R",
    output:
        "out/subsampling-experiment/summary-plot-historysizes.png",
        "out/subsampling-experiment/summary-plot-historysizes.svg",
        "out/subsampling-experiment/summary-plot-historysizes.rds"
    shell:
        """
        Rscript {input.rscript}
        """


rule plot_subsampling_r0_comparison:
    input:
        ss_orig_log,
        ss_mcmc_log_a,
        ss_mcmc_log_b,
        rscript="R/postprocessing-3-subsampling.R",
    output:
        "out/subsampling-experiment/summary-plot-r0.png",
        "out/subsampling-experiment/summary-plot-r0.svg",
        "out/subsampling-experiment/summary-plot-r0.rds"
    shell:
        """
        Rscript {input.rscript}
        """


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


rule run_subsample_mcmc:
    input:
       "out/subsampling-experiment/xml/{mcmc_name}.xml",
       beast_bin
    output:
        "out/{mcmc_name}.log"
    wildcard_constraints:
        mcmc_name = "timtam-2023-10-15-original|timtam-2023-10-15-subsample-0p66|timtam-2023-10-15-subsample-0p33"
    shell:
        """
        ./lib/beast/bin/beast -seed 1 -overwrite {input[0]}
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


rule make_ss_disaster_json:
    input:
        "config.xml",
        time_series_clean_csv,
        rscript = "R/preprocessing-3-subsampling.R"
    output:
        ss_disaster_json
    shell:
        """
        Rscript {input.rscript}
        """


rule plot_ss_disaster_png:
    input:
        "config.xml",
        ss_disaster_json,
        rscript = "R/preprocessing-4-subsampling.R"
    output:
        ss_disaster_png,
        ss_disaster_png.replace("png", "svg"),
        ss_disaster_png.replace("png", "rds"),
    shell:
        """
        Rscript {input.rscript}
        """


rule make_ss_mcmc_xml:
    input:
        "config.xml",
        ss_disaster_json,
        rscript = "R/preprocessing-5-subsampling.R"
    output:
        ss_orig_xml,
        ss_mcmc_xml_a,
        ss_mcmc_xml_b,
    shell:
        """
        Rscript {input.rscript}
        """
