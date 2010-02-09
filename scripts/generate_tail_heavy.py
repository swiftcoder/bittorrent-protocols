#! /usr/bin/env python 

import math, random

def generate(count, low, high, shape, func):
	result = [func(shape, i) for i in range(count)]
	
	L = min(result)
	H = max(result)
	
	for i in range(count):
		result[i] = (result[i] - L)/(H - L) * (high - low) + low
	
	return result

def generate_pareto(count, low, high, shape):
	return generate(count, low, high, shape, lambda x, i: random.paretovariate(x))

def generate_log_normal(count, low, high, shape):
	return generate(count, low, high, shape, lambda x, i: random.lognormvariate(1.0, x))

def generate_weibull(count, low, high, shape):
	return generate(count, low, high, shape, lambda x, i: random.weibullvariate(1.0, x))


def generate_exponential(count, low, high, shape):
	return generate(count, low, high, shape, lambda x, i: math.exp(i))

generators = {
	'pareto': generate_pareto,
	'log_normal': generate_log_normal,
	'weibull': generate_weibull,
	'exponential': generate_exponential,
}

def print_usage():
	print '\n\tusage: generate_tail_heavy.py [-h] [-count c] [-low l] [-high h] generator-name shape\n\n\t\twhere generator-name is one of:\n'
	for k, v in generators.iteritems():
		print '\t\t\t%s' % (k)
	print '\n\t\tand shape is a float in the range [0.0, 1.0].\n'

count = 100
low = 40.0
high = 12000.0

import sys

args = sys.argv[1:]
while len(args):
	if args[0] == '-count':
		count = int(args[1])
		args = args[2:]
	elif args[0] == '-low':
		low = float(args[1])
		args = args[2:]
	elif args[0] == '-high':
		high = float(args[1])
		args = args[2:]
	elif args[0] == '-h':
		print_usage()
		sys.exit()
	else:
		break

if len(args) < 2:
	print_usage()
	sys.exit()

which = str(args[0])
args = args[1:]

if not which in generators:
	print_usage()
	sys.exit()

#print 'generating: [%f:%f]' % (low, high)

#result = generate_pareto( count, low, high, float(args[0]) )
result = generators[which]( count, low, high, float(args[0]) )

result.sort()

for R in result:
	print R

