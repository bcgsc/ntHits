import argparse
import btllib
import os
import random
import subprocess


def parse_args():
    args = argparse.ArgumentParser()
    args.add_argument("-n", type=int, help="Number of sequences", default=300)
    args.add_argument("-k", type=int, help="k-mer length", default=8)
    return args.parse_args()


def generate_data(num_sequences):
    sequences = []
    for _ in range(num_sequences):
        length = random.randint(250, 1000)
        sequences.append(''.join(random.choices('ACGT', k=length)))
    return sequences


def write_fa(data):
    writer = btllib.SeqWriter('data.fa')
    for i in range(len(data)):
        writer.write(str(i), '', data[i])


def run_ntcard(k):
    subprocess.run([
        'ntcard', '-t', '48', '-k',
        str(k), '-c', '64', '-p', 'freq', 'data.fa'
    ])
    os.rename(f'freq_k{k}.hist', 'freq.hist')


def get_counts(k):
    subprocess.run([
        'dsk', '-file', 'data.fa', '-kmer-size',
        str(k), '-out', 'dsk.h5', '-abundance-min', '1'
    ])
    subprocess.run(['dsk2ascii', '-file', 'dsk.h5', '-out', 'counts.txt'])
    os.remove('dsk.h5')


def main():
    args = parse_args()
    data = generate_data(args.n)
    write_fa(data)
    run_ntcard(args.k)
    get_counts(args.k)


if __name__ == '__main__':
    main()