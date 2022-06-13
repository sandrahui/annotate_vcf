import sys
import json
import time
import vcf
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from urllib.error import HTTPError

# Script to annotate variants from a vcf file with:
# 1. Depth of sequence coverage at the site of variation.
# 2. Number of reads supporting the variant.
# 3. Percentage of reads supporting the variant versus those supporting
#    reference reads.
# 4. Using the VEP hgvs API, get the gene of the variant, type of variation
#    (substitution, insertion, CNV, etc.) and their effect (missense, silent,
#    intergenic, etc.). The API documentation is available here:
#    https://rest.ensembl.org/#VEP
# 5. The minor allele frequency of the variant if available.
# 6. Any additional annotations that you feel might be relevant.

# I included chrom/pos/ref/alt to identify variants, and kept the qual and
# filter columns from the original vcf file to make filtering easy. Because each
# variant can 1) contain multiple alternate alleles, and 2) affect multiple
# transcripts of the same gene, any given set of genomic coordinates may appear
# multiple times to allow for all combinations of alternate alleles and
# transcripts. However, this means that some information is repeated between
# lines, which is less space efficient. Missing values are denoted with "." (a
# period). Note this script uses the grch37 ensembl endpoint specifically
# (rather than the default grch38 endpoint) in order to match the reference of
# the given vcf file, which was aligned to /data/ref_genome/human_g1k_v37.fasta.

# Annotations (and column names) are:
# 1. total_seq_depth = Depth of sequence coverage at the site of variation.
# 2. num_alt_reads = Number of reads supporting the variant.
# 3. perc_alt_reads = Percentage of reads supporting the variant
#    perc_ref_reads = Percentage of reads supporting reference reads.
# 4. From the VEP hgvs API:
#    gene = gene of the variant, "." if missing
#    variant_type = type of variation (substitution, insertion, CNV, etc.)
#    variant_consequences = variant effect (missense, silent, intergenic, etc.),
#                           comma separated if multiple predicted effects, "." if missing
# 5. minor_allele_freq = the minor allele frequency of the variant if available, "." if missing
# 6. id = dbSNP rsID value, "." if missing
#    population_frequencies = comma separated list of population source and name
#                             mapped to allele frequency, "." if missing
#    transcript_id = Ensembl canonical transcript ID, "." if missing
#    variant_impact = rating given for compatibility with other variant annotation
#                     tools (HIGH, MODERATE, LOW, MODIFIER), "." if missing

# Expected input: vcf file
# Expected output: prints tab separated annotations to stdout, where columns are:
# chrom	pos	ref	alt	qual	filter	variant_type	total_seq_depth	num_alt_reads	perc_alt_reads	perc_ref_reads	minor_allele_freq	id	population_frequencies	gene	transcript_id	variant_impact	variant_consequences

# usage:
#  python annotate_vcf.py /path/to/input/vcf > /path/to/variant/annotations

# EnsemblRestClient code (with the exception of the get_variant_info() function)
# is from https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client
# (accessed 6/11/22)


class EnsemblRestClient(object):
    """ Code (__init__() and perform_rest_action functions()) from
    https://github.com/Ensembl/ensembl-rest/wiki/Example-Python-Client to use
    the Ensembl REST API to access data about variants while gracefully handling
    rate limiting
    """

    def __init__(self, server='http://rest.ensembl.org', reqs_per_sec=15):
        self.server = server
        self.reqs_per_sec = reqs_per_sec
        self.req_count = 0
        self.last_req = 0

    def perform_rest_action(self, endpoint, hdrs=None, params=None):
        if hdrs is None:
            hdrs = {}

        if 'Content-Type' not in hdrs:
            hdrs['Content-Type'] = 'application/json'

        if params:
            endpoint += '?' + urlencode(params)

        data = None

        # check if we need to rate limit ourselves
        if self.req_count >= self.reqs_per_sec:
            delta = time.time() - self.last_req
            if delta < 1:
                time.sleep(1 - delta)
            self.last_req = time.time()
            self.req_count = 0

        try:
            request = Request(self.server + endpoint, headers=hdrs)
            response = urlopen(request)
            content = response.read()
            if content:
                data = json.loads(content)
            self.req_count += 1

        except HTTPError as e:
            # check if we are being rate limited by the server
            if e.code == 429:
                if 'Retry-After' in e.headers:
                    retry = e.headers['Retry-After']
                    time.sleep(float(retry))
                    self.perform_rest_action(endpoint, hdrs, params)
            else:
                sys.stderr.write(
                    'Request failed for {0}: Status code: {1.code} Reason: {1.reason}\n'.format(endpoint, e))

        return data

    def get_variant_info(self, chrom, start, ref_allele, alt_allele):
        """ Fetches variant effect prediction information about a given variant

        Parameters
        ----------
        chrom : str
            chromosome number
        start : str
            starting coordinate of variant
        ref_allele : str
            reference allele of variant
        alt_allele : str
            alternate allele of variant to fetch information about

        Returns
        -------
        dict
            dictionary representation of the entire json returned from the ensembl query
        """
        end = start + len(ref_allele) - 1
        variant_info = self.perform_rest_action(
            endpoint=f"/vep/human/region/{chrom}:{start}:{end}/{alt_allele}", params={"variant_class": 1})[0]
        return variant_info


def extract_variant_info(variant_info, alt):
    """ Extracts relevant information from variant_info

    Parameters
    ----------
    variant_info : dict
        dictionary of json from ensembl, expected from get_variant_info()
    alt : str
        alternate allele of interest
    
    Returns
    -------
    str : variant_type
        type of variant/variant classification (SNP, indel, CNV, etc). See
        https://uswest.ensembl.org/info/genome/variation/prediction/classification.html
        for all possible values
    list : coloc_variants_info_list
        list of dicts, where each dict has minor allele frequence, ID, and
        population frequency information for colocated variants (ie, known
        variants at this location)
    list : transcript_info_list
        list of dicts, where each dict has gene, transcript ID, impact score,
        and variant consequences in each transcript. See
        https://uswest.ensembl.org/info/genome/variation/prediction/predicted_data.html
        for all possible variant consequences
    """
    variant_type = variant_info["variant_class"]
    allele_string = variant_info["allele_string"]

    # if have colocated variants, extract information about each one
    coloc_variants_info_list = []
    if("colocated_variants" in variant_info):
        for coloc_variants in variant_info["colocated_variants"]:
            # safety check: does the allele string we passed in match the
            # colocated variant? if not, then skip
            if(allele_string in coloc_variants["allele_string"]):
                # extract allele frequency and ID information, if they exist
                minor_allele_freq = "."
                if("minor_allele_freq" in coloc_variants):
                    minor_allele_freq = coloc_variants["minor_allele_freq"]
                pop_freqs_str = "."
                if("frequencies" in coloc_variants):
                    pop_freqs = []
                    for pop, freq in coloc_variants["frequencies"][str(alt)].items():
                        pop_freqs.append(":".join([pop, str(freq)]))
                    pop_freqs_str = ",".join(pop_freqs)
                id = "."
                if("id" in coloc_variants):
                    id = coloc_variants["id"],

                # collect into a dictionary
                coloc_variants_info_list.append({
                    "minor_allele_freq": minor_allele_freq,
                    "id": id,
                    "pop_freqs": pop_freqs_str
                })

    # if processed all colocated variants or didn't have any colocated variants,
    # and didn't find a match to the variant of interest, use a null list
    if(len(coloc_variants_info_list) == 0):
        coloc_variants_info_list.append({
            "minor_allele_freq": ".",
            "id": ".",
            "pop_freqs": "."
        })

    # if have transcripts, extract information about each one
    transcript_info_list = []
    if("transcript_consequences" in variant_info):
        for transcript in variant_info["transcript_consequences"]:
            # safety check, is this the allele of interest? if something went
            # wrong and we fetched info for a different variant at this site,
            # then skip
            if alt != transcript["variant_allele"]:
                continue
            
            # collect useful information into a dictionary
            transcript_info_list.append({
                "gene": transcript["gene_symbol"],
                "transcript_id": transcript["transcript_id"],
                "impact": transcript["impact"],
                "consequences": ",".join(transcript["consequence_terms"])
            })
    # if processed all transcript consequences or didn't have any transcript
    # consequences and didn't find a match to the variant of interest, use a
    # null list
    if(len(transcript_info_list) == 0):
        transcript_info_list.append({
            "gene": ".",
            "transcript_id": ".",
            "impact": ".",
            "consequences": "."
        })
    return variant_type, coloc_variants_info_list, transcript_info_list


if __name__ == "__main__":
    if(len(sys.argv) == 2):
        vcf_filename = sys.argv[1]
    else:
        print(f"usage: python {sys.argv[0]} </path/to/vcf>")
        sys.exit(0)

    # given vcf file was mapped to grch37, so use the grch37 specific ensembl endpoint
    client = EnsemblRestClient(server="http://grch37.rest.ensembl.org/")

    # print header
    print(f"chrom\tpos\tref\talt\tqual\tfilter\tvariant_type\ttotal_seq_depth\tnum_alt_reads\tperc_alt_reads\tperc_ref_reads\tminor_allele_freq\tid\tpopulation_frequencies\tgene\ttranscript_id\tvariant_impact\tvariant_consequences")

    # read in vcf file
    with(open(vcf_filename) as vcf_file):
        vcf_reader = vcf.Reader(vcf_file)

        # for each record in the vcf file
        for record in vcf_reader:
           # get basic variant identifying info
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            qual = record.QUAL
            filter = "PASS"
            if(len(record.FILTER) > 0):
                filter = ",".join(record.FILTER)

            # get number of reads for each allele
            # NR: Number of reads covering variant location in this sample (allele specific)
            # NV: Number of reads containing variant in this sample (allele specific)
            seq_depth_list = record.samples[0]["NR"]
            num_alt_reads_list = record.samples[0]["NV"]

            # each variant can have multiple alt alleles, so need to process each of them
            alt_list = record.ALT
            for alt, seq_depth, num_alt_reads in zip(alt_list, seq_depth_list, num_alt_reads_list):
                perc_alt_reads = int(num_alt_reads) / int(seq_depth)

                # fetch information about this variant from VEP rest api
                variant_info = client.get_variant_info(chrom, pos, ref, alt)

                # extract relevant information
                variant_type, coloc_variants_info_list, transcript_info_list = extract_variant_info(
                    variant_info, alt)

                # each position could have multiple colocated variants, so print each of them
                for coloc_variant in coloc_variants_info_list:
                    minor_allele_freq = coloc_variant["minor_allele_freq"]
                    id = coloc_variant["id"]
                    pop_freqs = coloc_variant["pop_freqs"]

                    # each alt allele can affect multiple transcripts, so print each of them
                    for transcript in transcript_info_list:
                        gene = transcript["gene"]
                        transcript_id = transcript["transcript_id"]
                        impact = transcript["impact"]
                        consequences = transcript["consequences"]
                        print(f"{chrom}\t{pos}\t{ref}\t{alt}\t{qual}\t{filter}\t{variant_type}\t{seq_depth}\t{num_alt_reads}\t{perc_alt_reads}\t{1-perc_alt_reads}\t{minor_allele_freq}\t{id}\t{pop_freqs}\t{gene}\t{transcript_id}\t{impact}\t{consequences}")
