import struct, time

def write_conan_string(f, S):
	bs = struct.pack("I", len(S))
	f.write(bs)
	f.write(S)
	f.write(bs)

def write_conan_main_header(f, B):
	S = "L=%f\nN=%u\nrank=%u\ncosmology=constraint\n" % (B.L, B.N, B.dim)
	write_conan_string(f, S)

def write_conan_sub_header(f, B, name, data):
	S = "name=%s\nL=%f\nN=%u\nrank=%u\nshape=%s\ntype=%s\n" \
		% (name, B.L, B.N, B.dim, B.shape, data.dtype)
	write_conan_string(f, S)

def write_conan_history(f):
	S = "%u:python constraint field code\n" % int(time.time())
	write_conan_string(f, S)

def write_conan_data(f, data):
	S = data.tostring()
	write_conan_string(f, S)

def write_conan(f, B, data):
	write_conan_main_header(f, B)
	write_conan_history(f)
	for k, d in data.iteritems():
		write_conan_sub_header(f, B, k, d)
		write_conan_data(f, d)

# vim:ts=4:sw=4:tw=80
