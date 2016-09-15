import subprocess
import os

def bayehem(k):
    """
    Bayehem function input k-mer val
    """
    if k % 2 == 0:
        k = k+1

    job_id = "k{}".format(k)

    assembly_path = evaluate(k)

    score = likelihood(assembly_path, k)

    return score


def evaluate(k):

    assembler_path = "./bayehem/assembler/rnaSPAdes-0.1.1-Linux/bin/rnaspades.py -o"

    data_path = "--pe1-1 bayehem/test_data/reads.left.fq \
                 --pe1-2 bayehem/test_data/reads.right.fq \
                 --pe2-1 bayehem/test_data/reads2.left.fq \
                 --pe2-2 bayehem/test_data/reads2.right.fq"
    output_dir = "bayehem/temp/{}".format(k[0])

    sh = " ".join([assembler_path, output_dir, "-k {}".format(k[0]), data_path])

    print("Evaluating k={}".format(k[0]))

    with open(os.devnull, 'w') as devnull:
        print(sh)
        if subprocess.call(sh, shell=True,
                           stdout=devnull, stderr=devnull) is not 0:
            raise RuntimeError("Assembler failed with param {}".format(k))

    return output_dir + "/contigs.fasta"


def likelihood(assembly_path, k):

    rsem_path = "./bayehem/rsem/detonate-1.11/rsem-eval/rsem-eval-calculate-score"
    data_path = "--paired-end bayehem/test_data/reads.left.fq,bayehem/test_data/reads2.left.fq \
                bayehem/test_data/reads.right.fq,bayehem/test_data/reads2.right.fq"
    out_name = "bayehem/temp/{}/score".format(k[0])
    length = "76"

    sh = " ".join([rsem_path, data_path, assembly_path, out_name, length])

    print("Calculating Likelihood")

    with open(os.devnull, 'w') as devnull:
        print(sh)
        if subprocess.call(sh, shell=True,
                           stdout=devnull, stderr=devnull) is not 0:
            raise RuntimeError("RSEM failed with param {}".format(k))

    with open(out_name+ ".score") as fh:
        score = float(fh.readline().strip().split("\t")[-1])

    return score

def main(job_id, params):
    print 'Anything printed here will end up in the output directory for job #%d' % job_id
    print params
    return bayehem(params['k'])#, job_id)

if __name__=='__main__':
    print(bayehem(13, "test"))
