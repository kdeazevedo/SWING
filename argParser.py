import argparse
import numpy as np

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
parser_run.add_argument('-n',type=int,default=1000,choices=range(1, 1000001),
        metavar="[1-1000000]", help='Number of sampling')
parser_run.add_argument('--dist',default='normal', choices = ["uniform","normal"],
        help='Sampling distribution')
parser_run.add_argument('-a','--angle',type=float,default=np.pi/24,
        help='Rotation max/min angles (in radian, ex np.pi/24)')
parser_run.add_argument('--seed',default=argparse.SUPPRESS,
        help='Seed for random sampling. Int or seed file\'s path')
parser_run.add_argument('--no-minimizer', action='store_true',
        help='Launch minimizer after sampling')

# Parser for sampling
parser_samples = subparsers.add_parser('samples',
        parents=[parent_parser],formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        conflict_handler='resolve')
parser_samples.add_argument('-n',type=int,default=1000,choices=range(1,1000001),
        metavar="[1-1000000]", help='Number of sampling')
parser_samples.add_argument('-c','--config', required=True, default=argparse.SUPPRESS,
        help='Config file for sampling')
parser_samples.add_argument('--dist',default='normal', choices = ["uniform","normal"],
        help='Sampling distribution')
parser_samples.add_argument('-a','--angle',type=float,default=np.pi/24,
        help='Rotation max/min angles (in radian, ex np.pi/24)')
parser_samples.add_argument('--seed',default=argparse.SUPPRESS,
        help='Seed for random sampling. Int or seed file\'s path')
parser_samples.add_argument('--no-minimizer', action='store_true',
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
