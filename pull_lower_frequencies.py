import pysam
import seaborn as sns
import matplotlib.pyplot as plt

def process_vcf(vcf_file):
    het_freqs = dict()
    vcf = pysam.VariantFile(vcf_file)
    for sample in vcf.header.samples:
        het_freqs[sample] = []
        new_vcf = pysam.VariantFile(vcf_file)
        for record in new_vcf:
            sample_data = record.samples[sample]
            genotype = sample_data["GT"]
            genotype_quality = sample_data.get("GQ", 0)
            read_depths = sample_data.get("AD", [])
            total_reads = sum(read_depths)

            # Only consider sites with GQ > 30 and depth > 20
            if genotype_quality:
                if genotype_quality > 30 and total_reads > 20 and genotype == (0, 1):
                    if total_reads > 0:
                        # Calculate the fraction for the less common allele
                        fewer_reads = min(read_depths)
                        fraction = fewer_reads / total_reads
                        het_freqs[sample].append(fraction)
    for sample in het_freqs:
        print(f"Qualifying hets in {sample}: {len(het_freqs[sample])}")
    return het_freqs

def plot_densities(heterozygous_freqs):
    samples = heterozygous_freqs.keys()
    plt.figure(figsize=(12, 10))
    for sample in samples:
        if sample != "F1_J_S392_L007":
            sample_name_parts = sample.split("_")
            sample_name = f"{sample_name_parts[0]} {sample_name_parts[1]}"
            sns.kdeplot(heterozygous_freqs[sample], label=sample_name, color="grey")
        else:
            sns.kdeplot(heterozygous_freqs[sample], label="Triploid", color='red')
    plt.xlabel("Frequency of less commonly observed allele", fontsize=16)
    plt.ylabel("Density", fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.legend(loc="upper left", fontsize=14)
    plt.savefig("Triturus_het_allelic_coverage_frequencies.pdf", bbox_inches='tight')


if __name__ == "__main__":
    het_freqs = process_vcf("F1s.biallelic.SNPs.q100.passOnly.vcf")
    plot_densities(het_freqs)