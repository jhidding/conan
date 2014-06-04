#!/usr/bin/python
#
# This python bit should be able to read the binary files and
# convert them to something python/gnuplot/ifrit ed.
# can understand
#

# a .conan file contains pairs of two blocks
# a block is delimited by two 4 byte integers giving the size of the block
# (Fortran style)
# each pair contains a header and data
# the header is '\n' delimited <key>=<value> pairs, of which some crucial
# values can be interpreted by a python 'eval' statement.
# the data entry folowing the first header gives the creation history of the file
# the following block pairs are header-data pairs where the header contains a
# 'name=' tag that gives information on the content of the data

# using these routines the header will be given as a dictionary with string
# keys, also the routine will try to reshape the data given the 'shape' tag in
# the header

import struct, sys, re, numpy, datetime, time
import numpy.random
import numpy as np

def read_block(f):
 	s = f.read(8)
 	if s == '':
 		raise Exception("EOF")
 
 	size  = struct.unpack("Q", s)[0]
 	block = f.read(size)
 	check = struct.unpack("Q", f.read(8))[0]
 
 	if check != size:
 		raise Exception("illigal block in file %s" % f.name)
 
 	return block
 
def read_header(f):
 	R = re.compile("([A-Za-z0-9_]*)=(.*)\n")
 	head = read_block(f).decode()
 	return dict(R.findall(head))
 
def read_history(f):
 	R = re.compile("([0-9]*):(.*)\n")
 	head = read_block(f).decode()
 	return dict(R.findall(head))

def read_array(f):
	head = read_header(f)
	print(head, file=sys.stderr)
	if 'name' in head:
		print("# %s" % head['name'], file=sys.stderr)
		
	if not 'dtype' in head:
		sys.stderr.write("# main header:\n")
		for k, i in head.items():
			sys.stderr.write("#\t%s = %s\n" % (k, i))
			
		sys.stderr.write("# history:\n")
		hist = read_history(f)
		times = sorted(map(int, hist.keys()))
		for t in times:
			sys.stderr.write("#\t%s: %s\n" % (datetime.datetime.fromtimestamp(t), hist[str(t)]))

		sys.stderr.flush()

		return head, hist

	else:
		print("# normal data block")
		dtype = head['dtype']
		raw_data = read_block(f)

		if dtype[0] == '[':
			dtype = np.dtype(eval(dtype), align=True)
			print(dtype, " to ", len(raw_data), file=sys.stderr)

		data = numpy.fromstring(raw_data, dtype)

	return (head, data)

def write_conan_string(f, S):
	bs = struct.pack("Q", len(S))
	f.write(bs)
	f.write(S)
	f.write(bs)

def log2(n):
    i = 0
    while n/2 >= 1:
        n /= 2
        i += 1
    return i
    
def write_conan_main_header(f, B):
	S = ("dim=%u\nN=%u\nsize=%f\ncosmology=Planck-nbody\nmbits=%u\n" \
        % (B.dim, B.N, B.L, log2(B.N))).encode()
	write_conan_string(f, S)

def write_conan_sub_header(f, B, name, data):
	S = ("name=%s\nL=%f\nN=%u\nrank=%u\nshape=%s\ntype=%s\n" \
		% (name, B.L, B.N, B.dim, B.shape, data.dtype)).encode()
	write_conan_string(f, S)

def write_conan_history(f):
	S = ("%u:python converted thing\n" \
        % int(time.time())).encode()
	write_conan_string(f, S)

def write_conan_data(f, data):
	S = data.tostring()
	write_conan_string(f, S)

def write_conan(f, B, data):
	write_conan_main_header(f, B)
	write_conan_history(f)
	for k, d in data.items():
		write_conan_sub_header(f, B, k, d)
		write_conan_data(f, d)


