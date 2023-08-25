import argparse
parser = argparse.ArgumentParser()

parser.add_argument("square", type=int,
                    help="display a square of a given number")
parser.add_argument("-v", "--verbosity", type=int, choices=[0, 1, 2],
                    help="increase output verbosity")
args = parser.parse_args()
answer = args.square**2
if args.v == 2:
    print "the square of {} equals {}".format(args.square, answer)
elif args.v == 1:
    print "{}^2 == {}".format(args.square, answer)
else:
    print answer
