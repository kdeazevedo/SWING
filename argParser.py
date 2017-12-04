import argparse

# Define common arguments
parent_parser = argparse.ArgumentParser(add_help=False)
parent_parser.add_argument('-o', default='out',
        help='The directory where program\'s output sotres')
parent_parser.add_argument('-rec', required=True, default=argparse.SUPPRESS,
        help='Recepter\'s file path')
parent_parser.add_argument('-lig', required=True, default=argparse.SUPPRESS,
        help='Ligand\'s file path')

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(dest='cmd')

# Parser for running from A to Z
parser_run = subparsers.add_parser('run',
    parents=[parent_parser],formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser_run.add_argument('-n',type=int,default=1000,help='Number of sampling')
parser_run.add_argument('-d',default='uniform', choices = ["uniform","normal"],
        help='Sampling distribution')
parser_run.add_argument('--minimizer', action='store_true',
        help='Launch minimizer after sampling')

# Parser for sampling
parser_samples = subparsers.add_parser('samples',
        parents=[parent_parser],formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        conflict_handler='resolve')
parser_samples.add_argument('-n',type=int,default=1000,help='Number of sampling')
parser_samples.add_argument('-c','--config', required=True, default=argparse.SUPPRESS,
        help='Config file for sampling')
parser_samples.add_argument('-d',default='uniform', choices = ["uniform","normal"],
        help='Sampling distribution')
parser_samples.add_argument('--minimizer', action='store_true',
        help='Launch minimizer after sampling')

# Parser for alignement
parser_align = subparsers.add_parser('align',
        parents=[parent_parser],formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_align.add_argument('-c','--config', required=True, default=argparse.SUPPRESS,
        help='Config file for alignement')

# Parser for download from InterEvol
parser_dl = subparsers.add_parser('download',
        parents=[parent_parser],formatter_class=argparse.ArgumentDefaultsHelpFormatter)

if __name__ == "__main__":
    p = parser.parse_args(['-rec','t','samples','-lig','u'])
    print(p)
