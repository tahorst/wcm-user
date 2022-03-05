#!/usr/bin/env python
"""
Run from wcEcoli directory on sherlock or locally to get correct paths to output and scripts.

TODO:
    - add option to copy output files to a central dir
    - remove dupes
"""

import argparse
import os
import shutil

from wholecell.utils import filepath as fp
from wholecell.utils import parallelization


SHERLOCK_PATH = '/home/travis/scratch/wcEcoli_'
COMPILED_PATH = '/home/travis/vb-share/2022-growth-paper/wcm-plots/'

# (fig, panel, sherlock, path, script, plot, label, options, out path, out labels)
ANALYSIS = [
    (2, 'a', True, 'out/20220121.200016__Amino_acid_combinations_in_media_without_regulation_or_charging/', 'Variant', 'doubling_time_histogram', '6+gen-var0,3-', '--generation-path-range 6 25 --variant-path 0 3', 'plotOut', ['_trimmed']),
    (2, 'b', True, 'out/20220116.130915__Amino_acid_combinations_in_media/', 'Variant', 'doubling_time_histogram', '6+gen-', '--generation-path-range 6 25', 'plotOut', ['_trimmed']),
    (2, 'c', True, 'out/20220116.130915__Amino_acid_combinations_in_media/', 'Variant', 'growth_trajectory', '', '--generation-path-range 6 25', '', []),
    (2, 'c', True, 'out/20220118.121155__Add_one_amino_acid_shift/', 'Variant', 'growth_trajectory', '', '--generation-path-range 2 16', '', []),
    (2, 'c', True, 'out/20220121.200028__Remove_one_amino_acid_shift/', 'Variant', 'growth_trajectory', '', '--generation-path-range 4 8', '', []),
    (2, 'c', True, 'out/20220123.144433__Conditions_without_regulation_or_charging/', 'Variant', 'growth_trajectory', '', '', '', []),
    (2, 'c', True, 'out/20220123.144515__Conditions_with_regulation/', 'Variant', 'growth_trajectory', '', '', '', []),
    (2, 'c', False, 'user/growth-paper/', 'growth_rp_plot.py', '', '', '', '', ['groups-combined-growth-rp']),
    (2, 'e', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation/', 'Parca', 'aa_synthesis_pathways', '', '', 'kb_plot_out', ['']),
    (2, 'f', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation/', 'Parca', 'amino_acid_uptake_rates', '', '', 'kb_plot_out', ['_clean']),
    (2, 'g', True, 'out/20220118.102135__Down_and_up_shifts_with_no_regulation/', 'Cohort', 'growth_time_series', '', '', 'timelines_000027/plotOut', ['_fig2']),
    (2, 'h', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation/', 'Cohort', 'growth_time_series', '', '', 'timelines_000027/plotOut', ['_fig2']),
    (3, 'a', True, 'out/20220120.184812__ppGpp_sensitivity_-_no_mechanistic_transport/', 'Variant', 'growth_trajectory', '', '--generation-path-range 2 8', '', []),
    (3, 'a', False, 'user/growth-paper/', 'growth_rp_plot.py', '', '', '', '', ['ppgpp-combined-growth-rp']),
    (3, 'bcde', True, 'out/20220120.184812__ppGpp_sensitivity_-_no_mechanistic_transport/', 'Variant', 'ppgpp_conc', '2+gen-minimal-', '--variant-path-range 0 10 --generation-path-range 2 8', 'plotOut', ['', '_output', '_capacity', '_excess']),
    (3, 'fg', True, 'out/20220120.184812__ppGpp_sensitivity_-_no_mechanistic_transport/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8', '', []),
    (3, 'f', True, 'out/20220117.215105__ppGpp_limitations_-_low_ppGpp/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8', '', []),
    (3, 'f', False, 'user/growth-paper/', 'ppgpp_growth.py', '', '', '', '', ['low-ppgpp']),
    (3, 'g', True, 'out/20220120.060837__ppGpp_limitations_-_high_ppGpp/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8', '', []),
    (3, 'g', True, 'out/20220304.172940__ppGpp_limitations_with_ribosomes_at_high_ppGpp/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8', '', []),
    (3, 'g', True, 'out/20220304.172940__ppGpp_limitations_with_ribosomes_at_high_ppGpp,_no_ppGpp_translation_inhibition/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8', '', []),
    (3, 'g', False, 'user/growth-paper/', 'ppgpp_growth.py', '', '', '', '', ['high-ppgpp']),
    (4, 'a', False, 'user/growth-paper/remove-aa-inhib/', 'conc_ki.py', '', '', '', '', ['side-by-side-bar']),
    (4, 'b', False, 'user/growth-paper/remove-aa-inhib/', 'conc_ki.py', '', '', '', '', ['aa-ki-prediction']),
    (4, 'c', True, 'out/20220119.081756__Remove_amino_acid_inhibition/', 'Cohort', 'growth_time_series', 'wt-', '-v0', 'remove_aa_inhibition_000000/plotOut', ['_fig4_single']),
    (4, 'c', True, 'out/20220119.081756__Remove_amino_acid_inhibition/', 'Cohort', 'growth_time_series', 'leuA-', '-v4', 'remove_aa_inhibition_000004/plotOut', ['_fig4_single']),
    (5, '', True, 'out/20220116.130915__Amino_acid_combinations_in_media/', 'Variant', 'growth_trajectory', '', '--generation-path-range 6 25', '', []),
    (5, '', False, 'user/growth-paper/', 'growth_rp_plot.py', '', '', '', '', ['shifts-combined-growth-rp']),
    (5, 'a', True, 'out/20220124.082636__Down_and_up_shifts_with_regulation_-_no_mechanistic_transport/', 'Variant', 'growth_trajectory', '', '', 'plotOut', ['_trimmed']),
    (5, 'a', True, 'out/20220124.082636__Down_and_up_shifts_with_regulation_-_no_mechanistic_transport/', 'Cohort', 'growth_time_series', '', '', 'timelines_000027/plotOut', ['_fig5']),
    (5, 'b', True, 'out/20220117.190915__Down_and_up_shifts_without_mechanistic_translation_supply/', 'Variant', 'growth_trajectory', '', '', 'plotOut', ['_trimmed']),
    (5, 'b', True, 'out/20220117.190915__Down_and_up_shifts_without_mechanistic_translation_supply/', 'Cohort', 'growth_time_series', '', '', 'timelines_000027/plotOut', ['_fig5']),
    (5, 'c', True, 'out/20220129.063633__Down_and_up_shifts_without_ppGpp_regulation_-_no_mechanistic_transport/', 'Variant', 'growth_trajectory', '', '', 'plotOut', ['_trimmed']),
    (5, 'c', True, 'out/20220129.063633__Down_and_up_shifts_without_ppGpp_regulation_-_no_mechanistic_transport/', 'Cohort', 'growth_time_series', '', '', 'timelines_000027/plotOut', ['_fig5']),
    (6, 'a', True, 'out/20220117.060726__Down_and_up_shifts_without_ppGpp_regulation/', 'Cohort', 'growth_time_series', '', '', 'timelines_000027/plotOut', ['_fig6']),
    (6, 'b', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation/', 'Cohort', 'growth_time_series', '', '', 'timelines_000027/plotOut', ['_fig6']),
    # (6, 'a', True, 'out/20220117.060726__Down_and_up_shifts_without_ppGpp_regulation_-_no_mechanistic_transport/', 'Cohort', 'growth_time_series', '', '', 'timelines_000027/plotOut', ['_fig6']),
    # (6, 'b', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation_-_no_mechanistic_transport/', 'Cohort', 'growth_time_series', '', '', 'timelines_000027/plotOut', ['_fig6']),
    ]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cpus', type=int, default=1)
    parser.add_argument('-f', '--fig', type=int)
    parser.add_argument('-p', '--panel')
    parser.add_argument('--sherlock', action='store_true')
    parser.add_argument('--local', action='store_true')
    parser.add_argument('--copy', action='store_true')
    return parser.parse_args()

def run_analysis(args):
    pool = parallelization.pool(num_processes=args.cpus)
    results = []
    for (fig, panel, sherlock, path, script, plot, label, options, *_) in ANALYSIS:
        if (args.fig == fig or args.fig is None) and (args.panel is None or args.panel in panel):
            if sherlock and args.sherlock:
                if label:
                    label = f'-o {label}'
                cmd = f'python runscripts/manual/analysis{script}.py {path} -p {plot} {label} {options}'
                print(cmd)
                results.append(pool.apply_async(fp.run_cmdline, (cmd,), kwds=dict(timeout=None)))
            if not sherlock and args.local:
                cmd = f'{path}{script}'
                print(cmd)
                results.append(pool.apply_async(fp.run_cmdline, (cmd,), kwds=dict(timeout=None)))
    pool.close()
    pool.join()

    # Check for errors
    for result in results:
        if not result.successful():
            result.get()

def copy_results():
    for (fig, panel, sherlock, path, _, plot, label, _, out_path, out_labels) in ANALYSIS:
        if (args.fig == fig or args.fig is None) and (args.panel is None or args.panel in panel):
            dest_dir = os.path.join(COMPILED_PATH, f'fig-{fig}')
            if not os.path.exists(dest_dir):
                print(f'*** Creating {dest_dir}')
                os.makedirs(dest_dir)

            for out_label in out_labels:
                if sherlock:
                    src = os.path.join(SHERLOCK_PATH + path, out_path, label + plot + out_label + '.pdf')
                else:
                    src = os.path.join(path, 'out', out_label + '.pdf')

                output_file = os.path.basename(src)
                dest = os.path.join(dest_dir, f'{panel}-{output_file}')

                if os.path.exists(src) and os.path.isfile(src):
                    print(f'Copying from {src} to {dest}')
                    shutil.copy2(src, dest)
                else:
                    print(f'*** Warning for {fig}{panel}: {src} is not a file')


if __name__ == '__main__':
    args = parse_args()

    if args.local or args.sherlock:
        run_analysis(args)

    if args.copy:
        copy_results()