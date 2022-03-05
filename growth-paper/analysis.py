#!/usr/bin/env python
"""
Run from wcEcoli directory on sherlock or locally to get correct paths to output and scripts.

TODO:
    - add option to copy output files to a central dir
    - remove dupes
"""

import argparse

from wholecell.utils import filepath as fp
from wholecell.utils import parallelization


analysis = [
    (2, 'a', True, 'out/20220121.200016__Amino_acid_combinations_in_media_without_regulation_or_charging/', 'Variant', 'doubling_time_histogram', '2a-6+gen-var0,3-', '--generation-path-range 6 25 --variant-path 0 3'),  # _trimmed
    (2, 'b', True, 'out/20220116.130915__Amino_acid_combinations_in_media/', 'Variant', 'doubling_time_histogram', '2b-6+gen-', '--generation-path-range 6 25'),  # _trimmed
    (2, 'c', True, 'out/20220116.130915__Amino_acid_combinations_in_media/', 'Variant', 'growth_trajectory', '', '--generation-path-range 6 25'),
    (2, 'c', True, 'out/20220118.121155__Add_one_amino_acid_shift/', 'Variant', 'growth_trajectory', '', '--generation-path-range 2 16'),
    (2, 'c', True, 'out/20220121.200028__Remove_one_amino_acid_shift/', 'Variant', 'growth_trajectory', '', '--generation-path-range 4 8'),
    (2, 'c', True, 'out/20220123.144433__Conditions_without_regulation_or_charging/', 'Variant', 'growth_trajectory', '', ''),
    (2, 'c', True, 'out/20220123.144515__Conditions_with_regulation/', 'Variant', 'growth_trajectory', '', ''),
    (2, 'c', False, 'user/growth-paper/', 'growth_rp_plot.py', '', '', ''),  # groups-
    (2, 'e', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation/', 'Parca', 'aa_synthesis_pathways', '2e-', ''),
    (2, 'f', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation/', 'Parca', 'amino_acid_uptake_rates', '2f-', ''),  # _clean
    (2, 'g', True, 'out/20220118.102135__Down_and_up_shifts_with_no_regulation/', 'Cohort', 'growth_time_series', '2g-', ''),  # _fig2
    (2, 'h', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation/', 'Cohort', 'growth_time_series', '2h-', ''),  # _fig2
    (3, 'a', True, 'out/20220120.184812__ppGpp_sensitivity_-_no_mechanistic_transport/', 'Variant', 'growth_trajectory', '', '--generation-path-range 2 8'),
    (3, 'a', False, 'user/growth-paper/', 'growth_rp_plot.py', '', '', ''),  # ppgpp-
    (3, 'bcde', True, 'out/20220120.184812__ppGpp_sensitivity_-_no_mechanistic_transport/', 'Variant', 'ppgpp_conc', '2+gen-minimal-', '--variant-path-range 0 10 --generation-path-range 2 8'),  # -, _output, _capacity, _excess
    (3, 'fg', True, 'out/20220120.184812__ppGpp_sensitivity_-_no_mechanistic_transport/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8'),
    (3, 'f', True, 'out/20220117.215105__ppGpp_limitations_-_low_ppGpp/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8'),
    (3, 'f', False, 'user/growth-paper/', 'ppgpp_growth.py', '', '', ''),  # low-
    (3, 'g', True, 'out/20220120.060837__ppGpp_limitations_-_high_ppGpp/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8'),
    (3, 'g', True, 'out/20220304.172940__ppGpp_limitations_with_ribosomes_at_high_ppGpp/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8'),
    (3, 'g', True, 'out/20220304.172940__ppGpp_limitations_with_ribosomes_at_high_ppGpp,_no_ppGpp_translation_inhibition/', 'Variant', 'growth_trajectory', '2+gen', '--generation-path-range 2 8'),
    (3, 'g', False, 'user/growth-paper/', 'ppgpp_growth.py', '', '', ''),  # high-
    (4, 'a', False, 'user/growth-paper/remove-aa-inhib/', 'conc_ki.py', '', '', ''),  # side-by-side-bar
    (4, 'b', False, 'user/growth-paper/remove-aa-inhib/', 'conc_ki.py', '', '', ''),  # aa-ki-prediction
    (4, 'c', True, 'out/20220119.081756__Remove_amino_acid_inhibition/', 'Cohort', 'growth_time_series', 'wt-', '-v0'),
    (4, 'c', True, 'out/20220119.081756__Remove_amino_acid_inhibition/', 'Cohort', 'growth_time_series', 'leuA-', '-v4'),
    (5, '', True, 'out/20220116.130915__Amino_acid_combinations_in_media/', 'Variant', 'growth_trajectory', '', '--generation-path-range 6 25'),
    (5, '', False, 'user/growth-paper/', 'growth_rp_plot.py', '', '', ''),  # shifts-
    (5, 'a', True, 'out/20220124.082636__Down_and_up_shifts_with_regulation_-_no_mechanistic_transport/', 'Variant', 'growth_trajectory', '5a-', ''),  # _trimmed
    (5, 'a', True, 'out/20220124.082636__Down_and_up_shifts_with_regulation_-_no_mechanistic_transport/', 'Cohort', 'growth_time_series', '5a-', ''),  # _fig5
    (5, 'b', True, 'out/20220117.190915__Down_and_up_shifts_without_mechanistic_translation_supply/', 'Variant', 'growth_trajectory', '5b-', ''),  # _trimmed
    (5, 'b', True, 'out/20220117.190915__Down_and_up_shifts_without_mechanistic_translation_supply/', 'Cohort', 'growth_time_series', '5b-', ''),  # _fig5
    (5, 'c', True, 'out/20220129.063633__Down_and_up_shifts_without_ppGpp_regulation_-_no_mechanistic_transport/', 'Variant', 'growth_trajectory', '5c-', ''),  # _trimmed
    (5, 'c', True, 'out/20220129.063633__Down_and_up_shifts_without_ppGpp_regulation_-_no_mechanistic_transport/', 'Cohort', 'growth_time_series', '5c-', ''),  # _fig5
    (6, 'a', True, 'out/20220117.060726__Down_and_up_shifts_without_ppGpp_regulation/', 'Cohort', 'growth_time_series', '6a-', ''),  # _fig6
    (6, 'b', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation/', 'Cohort', 'growth_time_series', '6b-', ''),  # _fig6
    # (6, 'a', True, 'out/20220117.060726__Down_and_up_shifts_without_ppGpp_regulation_-_no_mechanistic_transport/', 'Cohort', 'growth_time_series', '6a-', ''),  # _fig6
    # (6, 'b', True, 'out/20220116.124249__Down_and_up_shifts_with_regulation_-_no_mechanistic_transport/', 'Cohort', 'growth_time_series', '6b-', ''),  # _fig6
    ]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--cpus', type=int, default=1)
    parser.add_argument('-f', '--fig', type=int)
    parser.add_argument('-p', '--panel')
    parser.add_argument('--sherlock', action='store_true')
    parser.add_argument('--local', action='store_true')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()

    pool = parallelization.pool(num_processes=args.cpus)
    results = []
    for (fig, panel, sherlock, path, script, plot, label, options) in analysis:
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
