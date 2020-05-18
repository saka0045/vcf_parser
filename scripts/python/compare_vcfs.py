#!/usr/bin/python3

import argparse
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-b", "--baselineVcf", dest="baseline_vcf_path", required=True,
        help="Path to the baseline VCF to compare to"
    )
    parser.add_argument(
        "-c", "--compareVcf", dest="compare_vcf_path", required=True,
        help="Path to the VCF that is being compared to the baseline"
    )
    parser.add_argument(
        "-o", "--outDir", dest="output_directory", required=True,
        help="Path to the output directory to save result files"
    )
    parser.add_argument(
        "-a", "--alleleFrequencyDifference", dest="allele_frequency_difference", type=float, required=True,
        help="Allele frequency difference. Greater than this number will represent a match. "
             "Requires decimal between 0.00 and 1.00"
    )

    args = parser.parse_args()

    baseline_vcf_path = os.path.abspath(args.baseline_vcf_path)
    compare_vcf_path = os.path.abspath(args.compare_vcf_path)
    output_directory = os.path.abspath(args.output_directory)
    allele_frequency_difference = args.allele_frequency_difference

    baseline_vcf = open(baseline_vcf_path, "r")
    compare_vcf = open(compare_vcf_path, "r")

    # Parse the VCFs
    baseline_vcf_dict = {}
    parse_vcf(baseline_vcf, baseline_vcf_dict)
    compare_vcf_dict = {}
    parse_vcf(compare_vcf, compare_vcf_dict)

    # Compare the VCFs
    different_gt_dict, exact_match_dict, not_in_baseline_vcf_dict, not_in_compare_vcf_dict, \
    same_position_different_variant_dict = compare_vcfs(allele_frequency_difference, baseline_vcf_dict,
                                                        compare_vcf_dict)

    # Write out results
    exact_match_file = open(output_directory + "/exact_match_variants.csv", "w")
    different_gt_file = open(output_directory + "/match_with_different_gt.csv", "w")
    same_position_no_match_file = open(output_directory + "/same_position_no_match.csv", "w")
    not_in_compare_vcf_file = open(output_directory + "/not_in_compare_vcf.csv", "w")
    not_in_baseline_vcf_file = open(output_directory + "/not_in_baseline_vcf.csv", "w")

    write_results_found_in_both_vcfs(baseline_vcf_dict, compare_vcf_dict, exact_match_dict, exact_match_file)
    write_results_found_in_both_vcfs(baseline_vcf_dict, compare_vcf_dict, different_gt_dict, different_gt_file)
    write_results_found_in_both_vcfs(baseline_vcf_dict, compare_vcf_dict, same_position_different_variant_dict,
                                     same_position_no_match_file)
    write_results_only_found_in_one_vcf(baseline_vcf_dict, not_in_compare_vcf_dict, not_in_compare_vcf_file)
    write_results_only_found_in_one_vcf(compare_vcf_dict, not_in_baseline_vcf_dict, not_in_baseline_vcf_file)

    baseline_vcf.close()
    compare_vcf.close()
    exact_match_file.close()
    different_gt_file.close()
    same_position_no_match_file.close()
    not_in_compare_vcf_file.close()
    not_in_baseline_vcf_file.close()


def write_results_only_found_in_one_vcf(vcf_dict, variant_dict, result_file):
    """
    Write out results for variants that were only found in one of the VCFs
    :param vcf_dict:
    :param variant_dict:
    :param result_file:
    :return:
    """
    result_file.write("CHROM,POS,REF,ALT,GT,DP,AF\n")
    for chromosome, positions in variant_dict.items():
        for position in positions:
            ref = vcf_dict[chromosome][position]["ref"]
            alt = vcf_dict[chromosome][position]["alt"]
            gt = vcf_dict[chromosome][position]["GT"]
            dp = str(vcf_dict[chromosome][position]["DP"])
            af = str(vcf_dict[chromosome][position]["AF"])
            results = [ref, alt, gt, dp, af]
            result_file.write(chromosome + "," + position + "," + ",".join(results) + "\n")


def write_results_found_in_both_vcfs(baseline_vcf_dict, compare_vcf_dict, variant_dict, result_file):
    """
    Write out the results for variants that were found in both VCFs
    :param baseline_vcf_dict:
    :param compare_vcf_dict:
    :param variant_dict:
    :param result_file:
    :return:
    """
    result_file.write("CHROM,POS,REF,ALT,GT,DP,AF,,CHROM,POS,REF,ALT,GT,DP,AF\n")
    for chromosome, positions in variant_dict.items():
        for position in positions:
            baseline_ref = baseline_vcf_dict[chromosome][position]["ref"]
            baseline_alt = baseline_vcf_dict[chromosome][position]["alt"]
            baseline_gt = baseline_vcf_dict[chromosome][position]["GT"]
            baseline_dp = str(baseline_vcf_dict[chromosome][position]["DP"])
            baseline_af = str(baseline_vcf_dict[chromosome][position]["AF"])
            baseline_results = [baseline_ref, baseline_alt, baseline_gt, baseline_dp, baseline_af]
            compare_ref = compare_vcf_dict[chromosome][position]["ref"]
            compare_alt = compare_vcf_dict[chromosome][position]["alt"]
            compare_gt = compare_vcf_dict[chromosome][position]["GT"]
            compare_dp = str(compare_vcf_dict[chromosome][position]["DP"])
            compare_af = str(compare_vcf_dict[chromosome][position]["AF"])
            compare_results = [compare_ref, compare_alt, compare_gt, compare_dp, compare_af]
            result_file.write(chromosome + "," + position + "," + ",".join(baseline_results) + ",,")
            result_file.write(chromosome + "," + position + "," + ",".join(compare_results) + "\n")


def compare_vcfs(allele_frequency_difference, baseline_vcf_dict, compare_vcf_dict):
    exact_match_dict = {}
    different_gt_dict = {}
    same_position_different_variant_dict = {}
    not_in_compare_vcf_dict = {}
    not_in_baseline_vcf_dict = {}
    for chromosome in baseline_vcf_dict.keys():
        for position, baseline_contents in baseline_vcf_dict[chromosome].items():
            baseline_ref = baseline_contents["ref"]
            baseline_alt = baseline_contents["alt"]
            baseline_gt = baseline_contents["GT"]
            baseline_af = baseline_contents["AF"]
            # If the same position is found in the compare_vcf
            if position in compare_vcf_dict[chromosome].keys():
                compare_ref = compare_vcf_dict[chromosome][position]["ref"]
                compare_alt = compare_vcf_dict[chromosome][position]["alt"]
                compare_gt = compare_vcf_dict[chromosome][position]["GT"]
                compare_af = compare_vcf_dict[chromosome][position]["AF"]
                difference_af = abs(baseline_af - compare_af)
                # If ALT allele and GT is the same and the AF is within the specified amount,
                # add the position to exact_match_dict
                if baseline_alt == compare_alt and baseline_gt == compare_gt and \
                        difference_af <= allele_frequency_difference:
                    if chromosome not in exact_match_dict.keys():
                        exact_match_dict[chromosome] = [position]
                    else:
                        exact_match_dict[chromosome].append(position)
                # If ALT allele is the same and GT is different but the difference in AF is within 0.1
                elif baseline_alt == compare_alt and difference_af <= allele_frequency_difference and \
                        baseline_gt != compare_gt:
                    if chromosome not in different_gt_dict.keys():
                        different_gt_dict[chromosome] = [position]
                    else:
                        different_gt_dict[chromosome].append(position)
                # If the ALT alleles don't match or if the difference in AF is greater than 0.1
                else:
                    if chromosome not in same_position_different_variant_dict.keys():
                        same_position_different_variant_dict[chromosome] = [position]
                    else:
                        same_position_different_variant_dict[chromosome].append(position)
            # If the position is not found in compare VCF
            else:
                if chromosome not in not_in_compare_vcf_dict.keys():
                    not_in_compare_vcf_dict[chromosome] = [position]
                else:
                    not_in_compare_vcf_dict[chromosome].append(position)
    # If the position is not found in baseline vcf
    for chromosome in compare_vcf_dict.keys():
        for position in compare_vcf_dict[chromosome].keys():
            try:
                if position not in baseline_vcf_dict[chromosome].keys():
                    if chromosome not in not_in_baseline_vcf_dict.keys():
                        not_in_baseline_vcf_dict[chromosome] = [position]
                    else:
                        not_in_baseline_vcf_dict[chromosome].append(position)
            # If the chromosome is not found in baseline_vcf_dict, add it straight to the not_in_baseline_vcf_dict
            except KeyError:
                if chromosome not in not_in_baseline_vcf_dict.keys():
                    not_in_baseline_vcf_dict[chromosome] = [position]
                else:
                    not_in_baseline_vcf_dict[chromosome].append(position)
    return different_gt_dict, exact_match_dict, not_in_baseline_vcf_dict, not_in_compare_vcf_dict, same_position_different_variant_dict


def parse_vcf(vcf, vcf_dict):
    """
    Parses the vcf and stores the results in a dictionary. Each variant is categorized under chromosome and position
    :param vcf:
    :param vcf_dict:
    :return:
    """
    for line in vcf:
        # Skip the headers
        if line.startswith("#"):
            continue
        line = line.rstrip()
        line_items = line.split("\t")
        chromosome = line_items[0]
        position = line_items[1]
        ref = line_items[3]
        alt = line_items[4]
        format_key = line_items[8].split(":")
        format_value = line_items[9].split(":")
        genotype_index = format_key.index("GT")
        genotype = format_value[genotype_index]
        allele_depth_index = format_key.index("AD")
        allele_depth = format_value[allele_depth_index].split(",")
        ref_depth = allele_depth[0]
        alt_depth = allele_depth[1]
        total_depth = int(ref_depth) + int(alt_depth)
        allele_frequency = int(alt_depth) / int(total_depth)
        if chromosome not in vcf_dict.keys():
            vcf_dict[chromosome] = {}
        if position not in vcf_dict[chromosome].keys():
            vcf_dict[chromosome][position] = {}
            vcf_dict[chromosome][position]["ref"] = ref
            vcf_dict[chromosome][position]["alt"] = alt
            vcf_dict[chromosome][position]["GT"] = genotype
            vcf_dict[chromosome][position]["DP"] = total_depth
            vcf_dict[chromosome][position]["AF"] = round(allele_frequency, 2)
        elif position in vcf_dict[chromosome].keys():
            print("Variant found in duplicate position " + chromosome + position)
            print("This script will not work with duplicate positions, skipping this variant")
            continue


if __name__ == "__main__":
    main()
