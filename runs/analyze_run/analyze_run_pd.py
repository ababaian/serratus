#!/usr/bin/env python3
import sys
import seaborn as sns
import pandas as pd
from matplotlib import pyplot as plt
sns.set()

def read_accessions(fname):
    acc = pd.read_sql('acc', 'sqlite:///' + fname, index_col='acc_id')

    # Add columns from JSON to dataframe.
    acc['Run'] = [row.sra_run_info['Run'] for i, row in acc.iterrows()]
    for param in ('spots', 'bases', 'size_MB'):
        acc[param] = [int(row.sra_run_info[param]) for i, row in acc.iterrows()]

    # Calculate runtime, and convert from timedelta to seconds.
    dl_runtime = acc.split_end_time - acc.split_start_time
    merge_runtime = acc.merge_end_time - acc.merge_start_time

    acc['dl_runtime'] = [dt.total_seconds() for dt in dl_runtime]
    acc['merge_runtime'] = [dt.total_seconds() for dt in merge_runtime]

    return acc

def read_blocks(fname):
    # Read Blocks, and join with accession information
    blocks = pd.read_sql('blocks', 'sqlite:///' + fname, index_col='block_id')

    # Compute Block runtime
    align_runtime = blocks.align_end_time - blocks.align_start_time
    blocks['align_runtime'] = [dt.total_seconds() for dt in align_runtime]

    return blocks


def create_block_runtime_plots(blocks):
    f, axes = plt.subplots(1, 2, sharey=True)
    #axes[0].title('Block runtime grouped by Accession')
    sns.catplot(x='acc_id', y='align_runtime', data=blocks, ax=axes[0])

    #axes[1].title('Block runtime grouped by paired / unpaired')
    sns.catplot(x='contains_paired', y='align_runtime', data=blocks, ax=axes[1])

    #plt.savefig('grouped_block_times.svg'); plt.clf()
    plt.show()

def create_acc_runtime_plots(acc):
    param = 'size_MB'
    plt.title('Runtime per MB per worker')
    sns.regplot(param, 'dl_runtime', acc, label="DL")
    sns.regplot(param, 'merge_runtime', acc, label="Merge")
    sns.regplot(param, 'align_runtime', acc, label="Align")
    plt.legend(loc='best')
    #plt.savefig('accesion_times.svg'); plt.clf()
    plt.show()

    plt.title('Blocks per MB')
    sns.regplot(param, 'block_count', acc)
    plt.show()
    

def main():
    try:
        fname = sys.argv[1]
    except IndexError:
        print('usage: {} <filename>')
        sys.exit(1)

    acc = read_accessions(fname)
    blocks = read_blocks(fname)

    # Join accessions table into blocks
    blocks = blocks.join(acc, rsuffix='_acc', on='acc_id')

    # Add total block runtime to accessions
    acc['align_runtime'] = blocks.groupby(['acc_id'])['align_runtime'].sum()
    acc['block_count'] = blocks.groupby(['acc_id'])['align_runtime'].count()

    print('Total Input Size:     {} GB'.format(acc['size_MB'].sum() / 1000))
    print('Number of accessions: {}'.format(acc[acc['state'] == 'merge_done'].shape[0]))
    print('Number of blocks:     {}'.format(blocks[blocks['state'] == 'done'].shape[0]))

    create_block_runtime_plots(blocks)
    create_acc_runtime_plots(acc)

if __name__ == '__main__':
    main()


