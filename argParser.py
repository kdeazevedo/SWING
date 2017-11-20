import argparse

# Define common arguments
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('-o', default='out',
        help='The directory where program\'s output sotres')
parent_parser.add_argument('-rec', required=True, default=argparse.SUPPRESS,
        help='Recepter\'s file path')

parser = argparse.ArgumentParser(
        parents=[parent_parser],formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers()


# Parser for sampling
parser_samples = subparsers.add_parser('samples',
        parents=[parent_parser],formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        conflict_handler='resolve')
parser_samples.add_argument('-lig', nargs='+',required=True, default=argparse.SUPPRESS,
        help='Ligand\'s file path')
parser_samples.add_argument('-n',default=1000,help='Number of sampling')
parser_samples.add_argument('--minimizer',default=True,help='Lauch minimizer after sampling')

# Parser for alignement
parser_align = subparsers.add_parser('align',
        parents=[parent_parser],formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_align.add_argument('-lig', required=True, default=argparse.SUPPRESS,
        help='Ligand\'s file path')
parser_align.add_argument('-i','--interologs', nargs='+', required=True, default=argparse.SUPPRESS,
        help='Ligand\'s file path')

# Parser for download from InterEvol
parser_dl = subparsers.add_parser('download',
        parents=[parent_parser],formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_dl.add_argument('-lig', required=True, default=argparse.SUPPRESS,
        help='Ligand\'s file path')
