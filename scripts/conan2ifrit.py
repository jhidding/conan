#!/usr/bin/env python3

import thulsa, sys, struct
import numpy as np

def write_ifrit_string(f, S):
	bs = struct.pack("I", len(S))
	f.write(bs)
	f.write(S)
	f.write(bs)

def write_ifrit(f, header, data):
	N = len(data)
	L = eval(header['L'])
	write_ifrit_string(f, struct.pack("I", N))
	write_ifrit_string(f, struct.pack("6f", 0, 0, 0, L, L, L))
	write_ifrit_string(f, data[:,0].astype('float32').tostring())
	write_ifrit_string(f, data[:,1].astype('float32').tostring())
	write_ifrit_string(f, data[:,2].astype('float32').tostring())

fn_in  = sys.argv[1]
fn_out = ".".join(fn_in.split('.')[:-1]) + ".bin"

fi = open(fn_in, 'rb')
(header, history) = thulsa.read_array(fi)
(subhd, pos_data) = thulsa.read_array(fi)
fi.close()

fo = open(fn_out, 'wb')
write_ifrit(fo, header, pos_data)
fo.close()

