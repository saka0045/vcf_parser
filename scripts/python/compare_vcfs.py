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

    args = parser.parse_args()

    baseline_vcf_path = os.path.abspath(args.baseline_vcf_path)
    compare_vcf_path = os.path.abspath(args.compare_vcf_path)
    output_directory = os.path.abspath(args.output_directory)

    baseline_vcf = open(baseline_vcf_path, "r")
    compare_vcf = open(compare_vcf_path, "r")

    # Parse the VCFs
    baseline_vcf_dict = {}
    parse_vcf(baseline_vcf, baseline_vcf_dict)
    compare_vcf_dict = {}
    parse_vcf(compare_vcf, compare_vcf_dict)

    # Compare the VCFs
    exact_match_dict = {}
    different_gt_dict = {}
    same_position_different_variant_dict = {}
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
                # If ALT allele and GT is the same, add the position to exact_match_dict
                if baseline_alt == compare_alt and baseline_gt == compare_gt:
                    if chromosome not in exact_match_dict.keys():
                        exact_match_dict[chromosome] = [position]
                    else:
                        exact_match_dict[chromosome].append(position)
                # If ALT allele is the same and GT is different but the difference in AF is within 0.1
                elif baseline_alt == compare_alt and difference_af <= 0.1 and baseline_gt != compare_gt:
                    if chromosome not in different_gt_dict.keys():
                        different_gt_dict[chromosome] = [position]
                    else:
                        different_gt_dict[chromosome].append(position)
                else:
                    if chromosome not in same_position_different_variant_dict.keys():
                        same_position_different_variant_dict[chromosome] = [position]
                    else:
                        same_position_different_variant_dict[chromosome].append(position)

    print(exact_match_dict)
    print(different_gt_dict)
    print(same_position_different_variant_dict)

    baseline_vcf.close()
    compare_vcf.close()


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
        allele_frequency = int(alt_depth) / (int(alt_depth) + int(ref_depth))
        if chromosome not in vcf_dict.keys():
            vcf_dict[chromosome] = {}
        if position not in vcf_dict[chromosome].keys():
            vcf_dict[chromosome][position] = {}
            vcf_dict[chromosome][position]["ref"] = ref
            vcf_dict[chromosome][position]["alt"] = alt
            vcf_dict[chromosome][position]["GT"] = genotype
            vcf_dict[chromosome][position]["AF"] = round(allele_frequency, 2)
        elif position in vcf_dict[chromosome].keys():
            print("Variant found in duplicate position " + chromosome + position)
            print("This script will not work with duplicate positions, skipping this variant")
            continue


if __name__ == "__main__":
    main()
